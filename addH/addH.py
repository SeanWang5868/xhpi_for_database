import os
import gzip
import subprocess
import gemmi
import pandas as pd
from tqdm import tqdm
from datetime import datetime

# Directory path
today = datetime.now().strftime('%Y-%m-%d')
output_dir = os.path.join(".", today)
missing_monomers_file = os.path.join(output_dir, "missing_monomers.csv")

def is_gzip_file(filepath):
    """Check if a file is a gzip file."""
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def remove_residue_from_structure(structure, residue_name):
    """Remove residues with the specified name from the structure."""
    for model in structure:
        for chain in model:
            residues_to_remove = [i for i, res in enumerate(chain) if res.name == residue_name]
            for i in reversed(residues_to_remove):
                del chain[i]
    return structure

def process_gz_file(filepath):
    """Process a single .gz file to add hydrogen atoms and handle missing monomers."""
    try:
        with gzip.open(filepath, 'rb') as file:
            uncompressed_content = file.read().decode('utf-8')
        cif_block = gemmi.cif.read_string(uncompressed_content).sole_block()
        structure = gemmi.make_structure_from_block(cif_block)
    except Exception as e:
        print(f"Error processing CIF content for file {filepath}: {e}")
        return None

    temp_cif = os.path.join(output_dir, "temp.cif")
    temp_output_cif = os.path.join(output_dir, "temp_h.cif")

    try:
        structure.make_mmcif_document().write_file(temp_cif)
    except Exception as e:
        print(f"Error writing temporary CIF file for {filepath}: {e}")
        return None

    command = ['gemmi', 'h', '--water', temp_cif, temp_output_cif]
    max_retry = 8

    for _ in range(max_retry):
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"Successfully added hydrogen atoms for file {filepath}.")
            break
        else:
            error_msg = result.stderr.strip()
            if 'Monomer not in the library' in error_msg:
                monomer_name = error_msg.split('Monomer not in the library: ')[-1].split('.')[0]
                missing_monomers_df.loc[len(missing_monomers_df)] = [filepath, monomer_name]
                structure = remove_residue_from_structure(structure, monomer_name)
                try:
                    structure.make_mmcif_document().write_file(temp_cif)
                except Exception as e:
                    print(f"Error rewriting CIF after removing monomer {monomer_name} for file {filepath}: {e}")
                    return None
            else:
                print(f"Failed to add hydrogens for file {filepath}: {error_msg}")
                return None

    if os.path.exists(temp_output_cif):
        return temp_output_cif
    else:
        print(f"Failed to add hydrogens for file {filepath} after {max_retry} retries")
        return None

def save_modified_structure(temp_output_cif, original_filepath):
    """Save the modified structure file and compress it to a .gz file."""
    try:
        modified_structure = gemmi.read_structure(temp_output_cif)
        relative_path = os.path.relpath(original_filepath, output_dir)
        output_cif_gz = os.path.join(output_dir, relative_path)
        os.makedirs(os.path.dirname(output_cif_gz), exist_ok=True)
        with open(temp_output_cif, 'rt') as temp_output_file:
            with gzip.open(output_cif_gz, 'wt') as gz_output:
                gz_output.writelines(temp_output_file)
        print(f"Successfully saved the modified structure to {output_cif_gz}.")
    except Exception as e:
        print(f"Error writing gzipped CIF to file {output_cif_gz}: {e}")

def main():
    missing_monomers_df = pd.DataFrame(columns=["file", "monomer"])
    gz_files = [os.path.join(dirpath, filename)
                for dirpath, _, filenames in os.walk(output_dir)
                for filename in filenames if filename.endswith('.gz') and not filename.startswith('._')]

    with tqdm(total=len(gz_files), desc="Processing files") as pbar:
        for filepath in gz_files:
            if not is_gzip_file(filepath):
                pbar.update(1)
                continue
            
            temp_output_cif = process_gz_file(filepath)
            if temp_output_cif:
                save_modified_structure(temp_output_cif, filepath)
                os.remove(temp_output_cif)
            
            pbar.update(1)

    missing_monomers_df.to_csv(missing_monomers_file, index=False)
