import os
import gzip
import subprocess
import gemmi
import pandas as pd
from pathlib import Path
import logging
import tempfile
from typing import List, Tuple, Optional
from time import time
from tqdm import tqdm
import config

# Configure logging (minimal console output, detailed file output)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler('process_cif.log'),  # Detailed logs to file
        logging.StreamHandler()  # Only errors to console
    ]
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)  # Console shows only errors
file_handler = logging.FileHandler('process_cif.log')
file_handler.setLevel(logging.INFO)  # File logs INFO and above
logger.addHandler(file_handler)

# Define process_gz_file (from fixed Block 6)
def process_gz_file(filepath: Path, output_dir: Path, max_retry: int = 8) -> Tuple[Optional[Path], List[Tuple[str, str]]]:
    missing_monomers = []
    try:
        with gzip.open(filepath, 'rt', encoding='utf-8') as file:
            cif_content = file.read()
        cif_block = gemmi.cif.read_string(cif_content).sole_block()
        structure = gemmi.make_structure_from_block(cif_block)
    except (gzip.BadGzipFile, UnicodeDecodeError, RuntimeError) as e:
        logger.error(f"Failed to parse CIF content for {filepath}: {e}")
        return None, missing_monomers

    def remove_residue_from_structure(structure: gemmi.Structure, residue_name: str) -> gemmi.Structure:
        for model in structure:
            for chain in model:
                residues_to_remove = [i for i, res in enumerate(chain) if res.name == residue_name]
                for i in reversed(residues_to_remove):
                    del chain[i]
        return structure

    with tempfile.NamedTemporaryFile(mode='w+t', suffix='.cif', delete=False) as temp_cif, \
         tempfile.NamedTemporaryFile(mode='w+t', suffix='.cif', delete=False) as temp_output_cif:
        temp_cif_path = Path(temp_cif.name)
        temp_output_cif_path = Path(temp_output_cif.name)
        try:
            structure.make_mmcif_document().write_file(str(temp_cif_path))
            temp_cif.close()
            command = ['gemmi', 'h',  str(temp_cif_path), str(temp_output_cif_path)]
            for attempt in range(max_retry):
                result = subprocess.run(command, capture_output=True, text=True)
                if result.returncode == 0:
                    logger.info(f"Successfully added hydrogen atoms for {filepath}")
                    relative_path = filepath.relative_to(filepath.parents[1])
                    output_cif_gz = output_dir / relative_path.with_suffix('.gz')
                    output_cif_gz.parent.mkdir(parents=True, exist_ok=True)
                    with open(temp_output_cif_path, 'rt', encoding='utf-8') as temp_output_file:
                        with gzip.open(output_cif_gz, 'wt', encoding='utf-8') as gz_output:
                            gz_output.writelines(temp_output_file)
                    logger.info(f"Saved modified structure to {output_cif_gz}")
                    return output_cif_gz, missing_monomers
                else:
                    error_msg = result.stderr.strip()
                    if 'Monomer not in the library' in error_msg:
                        monomer_name = error_msg.split('Monomer not in the library: ')[-1].split('.')[0]
                        missing_monomers.append((str(filepath), monomer_name))
                        logger.warning(f"Monomer {monomer_name} not found for {filepath}, removing and retrying")
                        structure = remove_residue_from_structure(structure, monomer_name)
                        try:
                            structure.make_mmcif_document().write_file(str(temp_cif_path))
                        except RuntimeError as e:
                            logger.error(f"Failed to rewrite CIF after removing {monomer_name} for {filepath}: {e}")
                            return None, missing_monomers
                    else:
                        logger.error(f"Failed to add hydrogens for {filepath}: {error_msg}")
                        return None, missing_monomers
            logger.error(f"Failed to add hydrogens for {filepath} after {max_retry} retries")
            return None, missing_monomers
        except (subprocess.SubprocessError, RuntimeError, IOError) as e:
            logger.error(f"Error processing {filepath}: {e}")
            return None, missing_monomers
        finally:
            try:
                temp_cif_path.unlink(missing_ok=True)
                temp_output_cif_path.unlink(missing_ok=True)
            except OSError as e:
                logger.warning(f"Failed to clean up temporary files for {filepath}: {e}")

# Set up variables
input_dir = Path('/media/bql506/Sean/test')
output_dir = Path('/media/bql506/Sean/test_addH')
max_retry = 3

if not input_dir.is_dir():
    logger.error(f"Input directory {input_dir} does not exist or is not a directory")
    raise SystemExit(1)
output_dir.mkdir(parents=True, exist_ok=True)

gz_files = [f for f in input_dir.rglob('*.gz') if not f.name.startswith('._')]
logger.info(f"Found {len(gz_files)} .gz files to process")
print(f"Processing {len(gz_files)} files")

# Process files with tqdm progress bar
missing_monomers = []
processed_count = 0
failed_count = 0
total_files = len(gz_files)
start_time = time()
avg_time_per_file = 0

pbar = tqdm(total=total_files, desc="Processing files", unit="file")

for filepath in gz_files:
    file_start_time = time()
    try:
        output_cif, monomers = process_gz_file(filepath, output_dir, max_retry)
        missing_monomers.extend(monomers)
        processed_count += 1
        logger.info(f"Processed {filepath}, output: {output_cif}")
    except Exception as e:
        failed_count += 1
        logger.error(f"Unexpected error processing {filepath}: {e}")

    # Update average time
    elapsed = time() - file_start_time
    avg_time_per_file = ((avg_time_per_file * (processed_count + failed_count - 1)) + elapsed) / (processed_count + failed_count)

    # Update tqdm progress bar
    pbar.update(1)
    pbar.set_postfix({"Processed": processed_count, "Failed": failed_count})

pbar.close()

# Clear the updating line and print final summary
print(f"Completed: {processed_count}/{len(gz_files)} files processed, {failed_count} failed")

# Save missing monomers to a local CSV file
missing_monomers_file = output_dir / "missing_monomers.csv"
if missing_monomers:
    missing_monomers_df = pd.DataFrame(missing_monomers, columns=["file", "monomer"])
    missing_monomers_df.to_csv(missing_monomers_file, index=False)
    logger.info(f"Saved missing monomers to {missing_monomers_file}")
else:
    logger.info("No missing monomers found")

# Load CSV to compute summary (to save memory as per request)
if missing_monomers_file.exists():
    df = pd.read_csv(missing_monomers_file)
    total_missing_residues = len(df)
    total_missing_pdbs = df['file'].nunique()
    print(f"Total missing PDBs: {total_missing_pdbs}")
    print(f"Total missing residues: {total_missing_residues}")
else:
    print("Total missing PDBs: 0")
    print("Total missing residues: 0")