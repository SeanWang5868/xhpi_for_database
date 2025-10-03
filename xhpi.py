import os
import sys
import pandas as pd
import gzip
import gemmi
from tqdm import tqdm
from datetime import datetime
import wget
from calculate import printer,detect_plevin
from addH import addH

ring_atoms_dict = {
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'],
    'PHE': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'],
    'HIS': ['CE1', 'ND1', 'NE2', 'CG', 'CD2']}
trp_A_dict = {
    'TRP': ['CD1', 'CD2', 'NE1', 'CG', 'CE2']}

printer.print_xhpi()

# Create a folder with today's date as the name
today = datetime.now().strftime('%Y-%m-%d')
output_dir = os.path.join("../output_dir", today)
os.makedirs(output_dir, exist_ok=True)

print(f"Output directory created: {output_dir}")

# Receive 4-letter PDB names for download
pdb_names = []
print("Enter 4-letter PDB names separated by commas (e.g., ABCD,EFGH,IJKL):")
input_str = input().strip().upper()

# Process the input string, remove spaces, and split by comma
pdb_names = [name for name in input_str.split(',') if len(name) == 4]

if len(pdb_names) == 0:
    print("No valid 4-letter PDB names provided. Exiting.")
    sys.exit(1)

print(f"PDB to be processed: {', '.join(pdb_names)}")

# Download PDB files
for pdb in tqdm(pdb_names, desc="Downloading PDB files"):
    url = f"https://files.rcsb.org/download/{pdb}.cif.gz"
    output_path = os.path.join(output_dir, f"{pdb}.cif.gz")
    
    if os.path.exists(output_path):
        print(f"\n{pdb} already exists in {output_path}. Skipping download.")
        continue
    try:
        wget.download(url, output_path)
        print(f"\nDownloaded and saved {pdb} to {output_path}")
    except Exception as e:
        print(f"Error downloading {pdb}: {e}")

# Process the downloaded files and calculate XHPI
gz_files = [os.path.join(output_dir, f"{pdb_name}.cif.gz") for pdb_name in pdb_names]

addH.main()

result = []

with tqdm(total=len(gz_files), desc="Processing files") as pbar:
    for filepath in gz_files:
        try:
            pdb_name = os.path.basename(filepath).replace('.cif.gz', '')

            # Read the gzipped CIF file
            with gzip.open(filepath, 'rb') as file:
                uncompressed_content = file.read().decode('utf-8')
                cif = gemmi.cif.read_string(uncompressed_content).sole_block()
                structure = gemmi.make_structure_from_block(cif)
                resolution = structure.resolution
                model = structure[0]

                for chain in model:
                    for residue in chain:
                        
                        if residue.name in ring_atoms_dict:
                            found_interactions = detect_plevin.detect_plevin(
                                pdb_name, resolution, model, chain, structure, residue)
                            result.extend(found_interactions)

                        if residue.name in trp_A_dict:
                            found_interactions = detect_plevin.detect_plevin(
                                pdb_name, resolution, model, chain, structure, residue)
                            result.extend(found_interactions)

        except Exception as e:
            print(f"Error processing file {filepath}: {e}")

        finally:
            pbar.update(1)

print(f"Detected {len(result)} XH-π interactions in total.")

# Save the calculation results to a CSV file
if result: # Only save if results were found
    df = pd.DataFrame(result)
    output_path = os.path.join(output_dir, 'xhpi_output.csv')
    df.to_csv(output_path, index=False)
    print(f'The result has been saved to {output_path}')
else:
    print("No interactions were found to save.")


