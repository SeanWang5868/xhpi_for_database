import os
import sys
import pandas as pd
import gzip
import gemmi
import config
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

input_directory = "/Users/bql506/Desktop/xhpi_for_database/pdb_addedH_u7"
# input_directory = "/media/bql506/Sean/test_addH"

gz_files = []
for dirpath, _, filenames in os.walk(input_directory):
    for filename in filenames:
        if filename.endswith(".cif.gz"):
            # 将文件的完整路径添加到列表中
            gz_files.append(os.path.join(dirpath, filename))
gz_files



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
                            found_interactions = detect_plevin.detect_plevin(pdb_name, resolution, model, chain, structure, residue, ring_atoms_dict)
                            result.extend(found_interactions)     

                        if residue.name in trp_A_dict:
                            found_interactions = detect_plevin.detect_plevin(pdb_name, resolution, model, chain, structure, residue, trp_A_dict)
                            result.extend(found_interactions)

        except Exception as e:
            print(f"Error processing file {filepath}: {e}")

        finally:
            pbar.update(1)

xhpi_true_count = sum(1 for interaction in result if interaction.get('xhpi') == 1)
print(f"Found {xhpi_true_count} interactions meeting the XH-π criteria")

# Save the calculation results to a CSV file
if result: # Only save if results were found
    df = pd.DataFrame(result)
    output_path = os.path.join(input_directory, 'xhpi_output.csv')
    df.to_csv(output_path, index=False)
    print(f'The result has been saved to {output_path}')
else:
    print("No interactions were found to save.")