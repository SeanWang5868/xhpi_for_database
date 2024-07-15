import os
import sys
from tqdm import tqdm
from datetime import datetime
import wget
import os
import gzip
import subprocess
import gemmi
import pandas as pd
from calculate import get_info, calculator, append_result, rename, printer, Plevin, Plevin_append_result

printer.print_xhpi()

# create a folder
today = datetime.now().strftime('%Y-%m-%d')
output_dir = os.path.join(".", today)
os.makedirs(output_dir, exist_ok=True)

print(f"Output directory created: {output_dir}")

# 4 charactors PDB name
pdb_names = []
print("Enter 4-letter PDB names separated by commas (e.g., 12as,5FJJ,1Gn0):")
input_str = input("Enter PDB names: ").strip().upper()
pdb_names = [name for name in input_str.split(',') if len(name) == 4]

if len(pdb_names) == 0:
    print("No valid 4-letter PDB names provided. Exiting.")
    sys.exit(1)

print(f"PDB names to be processed: {', '.join(pdb_names)}")

# save PDB files to folder
for pdb in tqdm(pdb_names, desc="Downloading PDB files"):
    url = f"https://files.rcsb.org/download/{pdb}.cif.gz"
    output_path = os.path.join(output_dir, f"{pdb}.cif.gz")
    try:
        wget.download(url, output_path)
        print(f"\nDownloaded and saved {pdb} to {output_path}")
    except Exception as e:
        print(f"Error downloading {pdb}: {e}")

# choose the method: Hudson or Plevin
method = input("Enter method (Hudson or Plevin): ").strip().lower()

if method not in ['hudson', 'plevin']:
    print("Invalid method. Please enter 'Hudson' or 'Plevin'.")
    sys.exit(1)

gz_files = [os.path.join(output_dir, f"{pdb_name}.cif.gz") for pdb_name in pdb_names]

# set the path of Monomer library
os.environ['CLIBD_MON'] = '/Volumes/Sean/pdb_mirror/monomer/monomers_17June'

root_dir = today
# root_dir = "/Users/seanwang/Library/CloudStorage/GoogleDrive-bql506@york.ac.uk/My Drive/AlphaFold/pdb_mirror/data/structures/divided"
# output_dir = "/Volumes/Sean/pdb_mirror/Add_H/added_H_Test"
missing_monomers_file = os.path.join(output_dir, "missing_monomers.csv")

def is_gzip_file(filepath):
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

# the df to save the missing monomers
missing_monomers_df = pd.DataFrame(columns=["file", "monomer"])

def remove_residue_from_structure(structure, residue_name):
    for model in structure:
        for chain in model:
            residues_to_remove = [i for i, res in enumerate(chain) if res.name == residue_name]
            for i in reversed(residues_to_remove):
                del chain[i]  # use __delitem__ to delete the missing residue
    return structure

# find all gz. files
gz_files = []
for dirpath, _, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.endswith('.gz') and not filename.startswith('._'):
            gz_files.append(os.path.join(dirpath, filename))

with tqdm(total=len(gz_files), desc="Processing files") as pbar:
    for filepath in gz_files:
        if not is_gzip_file(filepath):
            pbar.update(1)
            continue
        
        try:
            with gzip.open(filepath, 'rb') as file:
                uncompressed_content = file.read().decode('utf-8')
            cif_block = gemmi.cif.read_string(uncompressed_content).sole_block()
            structure = gemmi.make_structure_from_block(cif_block)
        except Exception as e:
            print(f"Error processing CIF content for file {filepath}: {e}")
            pbar.update(1)
            continue

        temp_cif = os.path.join(output_dir, "temp.cif")
        structure.make_mmcif_document().write_file(temp_cif)
        
        temp_output_cif = os.path.join(output_dir, "temp_h.cif")
        command = ['gemmi', 'h', temp_cif, temp_output_cif]
        max_retry = 8  # max trying times 
        retry_count = 0
        success = False
        while retry_count < max_retry:
            retry_count += 1
            result = subprocess.run(command, capture_output=True, text=True)
            
            if result.returncode == 0:
                success = True
                break
            else:
                error_msg = result.stderr.strip()
                if 'Monomer not in the library' in error_msg:
                    monomer_name = error_msg.split('Monomer not in the library: ')[-1].split('.')[0]
                    new_row = pd.DataFrame({"file": [filepath], "monomer": [monomer_name]})
                    missing_monomers_df = pd.concat([missing_monomers_df, new_row], ignore_index=True)
                    try:
                        structure = remove_residue_from_structure(structure, monomer_name)
                        structure.make_mmcif_document().write_file(temp_cif)
                    except Exception as e:
                        print(f"Error processing structure after removing monomer {monomer_name} for file {filepath}: {e}")
                        pbar.update(1)
                        continue
                else:
                    print(f"Failed to add hydrogens for file {filepath}")
                    pbar.update(1)
                    break

        if success:
            try:
                modified_structure = gemmi.read_structure(temp_output_cif)
                cif_string = modified_structure.make_mmcif_document().as_string()
                relative_path = os.path.relpath(filepath, root_dir)
                output_cif_gz = os.path.join(output_dir, relative_path)
                output_cif_gz_dir = os.path.dirname(output_cif_gz)
                if not os.path.exists(output_cif_gz_dir):
                    os.makedirs(output_cif_gz_dir)
                with open(temp_output_cif, 'rt') as temp_output_file:
                    with gzip.open(output_cif_gz, 'wt') as gz_output:
                        gz_output.writelines(temp_output_file)
            except Exception as e:
                print(f"Error writing gzipped CIF to file {output_cif_gz}: {e}")
                pbar.update(1)
        else:
            print(f"Failed to add hydrogens for file {filepath} after {max_retry} retries")
            pbar.update(1)

        os.remove(temp_cif)
        os.remove(temp_output_cif)
        pbar.update(1)
        missing_monomers_df.to_csv(missing_monomers_file, index=False)

print("Successfully completed hydrogen adding.")

ring_atoms_dict = {
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'],
    'PHE': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CG'],
    'HIS': ['CE1', 'ND1', 'NE2', 'CG', 'CD2']}
trp_A_dict = {
    'TRP': ['CD1', 'CD2', 'NE1', 'CG', 'CE2']}

result = []

with tqdm(total=len(gz_files), desc="Processing files") as pbar:
    for filepath in gz_files:
        try:
            pdb_name = os.path.basename(filepath).replace('.cif.gz', '')  # get pdb name

            # read all .cif.gz files
            with gzip.open(filepath, 'rb') as file:
                uncompressed_content = file.read().decode('utf-8')
                cif = gemmi.cif.read_string(uncompressed_content).sole_block()
                structure = gemmi.make_structure_from_block(cif)
                resolution = structure.resolution

                neighbor_search = gemmi.NeighborSearch(structure[0], structure.cell, 6)
                neighbor_search.populate(include_h=True)

                model = structure[0]
                for chain in model:
                    for residue in chain:
                        if residue.name in ring_atoms_dict:
                            pi_atoms = [atom for atom in residue if atom.name in ring_atoms_dict[residue.name]]

                            if len(pi_atoms) == len(ring_atoms_dict[residue.name]):
                                pi_center, pi_center_array, pi_atoms_array, normal_vector, pi_b_factor_mean = get_info.piinfo(pi_atoms)
                                # find X atoms
                                alt_pi = pi_atoms[0].altloc
                                X_atoms = neighbor_search.find_atoms(pi_center, alt=alt_pi, min_dist=0.0, radius=6)

                                if method == 'hudson':
                                    # Hudson System
                                    for X_atom in X_atoms:
                                        X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
                                        mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2

                                        if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                                            X_to_pi_center_distance = calculator.distance(X_pos_array, pi_center_array)

                                            if X_to_pi_center_distance < 6:
                                                alt_H = X_atom.altloc
                                                H_atoms = neighbor_search.find_atoms(X_cra.atom.pos, alt=alt_H, min_dist=0.0, radius=1.3)
                                                for H_atom in H_atoms:
                                                    H_cra, H_chain_name, H_residue_num, H_residue_name, H_element_name, H_atom_name, H_b_iso, H_pos, H_pos_array = get_info.atominfo(H_atom, model)
                                                    
                                                    if H_cra.atom.element == gemmi.Element('H'):
                                                        angle_in_degrees = calculator.theta_angle(pi_center_array, X_pos_array, H_pos_array, normal_vector)
                                                        
                                                        if angle_in_degrees is not None and angle_in_degrees <= 40.0:
                                                            projection_distance = calculator.projection_distance(normal_vector, pi_center_array, X_pos_array)

                                                            if (residue.name == 'HIS' and projection_distance <= 1.6) or (residue.name in ['TRP', 'TYR', 'PHE'] and projection_distance <= 2.0):
                                                                pi_residue_name = rename.pi_rename(residue)

                                                                append_result.append_result(result, pdb_name, resolution, chain, residue, mean_b_factor,
                                                                    pi_residue_name, pi_center_array, normal_vector, pi_b_factor_mean, 
                                                                    X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_pos_array, X_b_iso,
                                                                    H_atom_name, H_pos_array, X_to_pi_center_distance, angle_in_degrees, projection_distance)

                                elif method == 'plevin':
                                    # Plevin System
                                    for X_atom in X_atoms:
                                        X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
                                        mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2

                                        if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                                            X_to_pi_center_distance = Plevin.distance(X_pos_array, pi_center_array)

                                            if X_to_pi_center_distance < 6:
                                                XPCN_angle = Plevin.XPCN_angle(X_pos_array, pi_center_array, normal_vector)

                                                if XPCN_angle < 25.0:
                                                    alt_H = X_atom.altloc
                                                    H_atoms = neighbor_search.find_atoms(X_cra.atom.pos, alt=alt_H, min_dist=0.0, radius=1.3)
                                                    for H_atom in H_atoms:
                                                        H_cra, H_chain_name, H_residue_num, H_residue_name, H_element_name, H_atom_name, H_b_iso, H_pos, H_pos_array = get_info.atominfo(H_atom, model)
                                                        if H_cra.atom.element == gemmi.Element('H'):
                                                            XH_picenterAngle = Plevin.XH_picenterAngle(pi_center_array, X_pos_array, H_pos_array) 
                                                            
                                                            if XH_picenterAngle is not None and XH_picenterAngle > 120.0:                        
                                                                pi_residue_name = 'TRP_A'

                                                                Plevin_append_result.append_result(result, pdb_name, resolution, chain, residue, mean_b_factor,
                                                                    pi_residue_name, pi_center_array, normal_vector, pi_b_factor_mean, 
                                                                    X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_pos_array, X_b_iso,
                                                                    H_atom_name, H_pos_array, X_to_pi_center_distance, XH_picenterAngle, XPCN_angle)

                        if residue.name in trp_A_dict:
                            pi_atoms = [atom for atom in residue if atom.name in trp_A_dict[residue.name]]

                            if len(pi_atoms) == len(trp_A_dict[residue.name]):
                                pi_center, pi_center_array, pi_atoms_array, normal_vector, pi_b_factor_mean = get_info.piinfo(pi_atoms)
                                # Find X atom
                                alt_pi = pi_atoms[0].altloc
                                X_atoms = neighbor_search.find_atoms(pi_center, alt=alt_pi, min_dist=0.0, radius=6)

                                if method == 'hudson':
                                    for X_atom in X_atoms:
                                        X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
                                        mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2

                                        if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                                            X_to_pi_center_distance = calculator.distance(X_pos_array, pi_center_array)

                                            if X_to_pi_center_distance < 6:
                                                alt_H = X_atom.altloc
                                                H_atoms = neighbor_search.find_atoms(X_cra.atom.pos, alt=alt_H, min_dist=0.0, radius=1.3)
                                                for H_atom in H_atoms:
                                                    H_cra, H_chain_name, H_residue_num, H_residue_name, H_element_name, H_atom_name, H_b_iso, H_pos, H_pos_array = get_info.atominfo(H_atom, model)
                                                    
                                                    if H_cra.atom.element == gemmi.Element('H'):
                                                        angle_in_degrees = calculator.theta_angle(pi_center_array, X_pos_array, H_pos_array, normal_vector)
                                                        
                                                        if angle_in_degrees is not None and angle_in_degrees <= 40.0:
                                                            projection_distance = calculator.projection_distance(normal_vector, pi_center_array, X_pos_array)

                                                            if (residue.name == 'HIS' and projection_distance <= 1.6) or (residue.name in ['TRP', 'TYR', 'PHE'] and projection_distance <= 2.0):
                                                                pi_residue_name = 'TRP_A'

                                                                append_result.append_result(result, pdb_name, resolution, chain, residue, mean_b_factor,
                                                                    pi_residue_name, pi_center_array, normal_vector, pi_b_factor_mean, 
                                                                    X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_pos_array, X_b_iso,
                                                                    H_atom_name, H_pos_array, X_to_pi_center_distance, angle_in_degrees, projection_distance)

                                elif method == 'plevin':
                                    for X_atom in X_atoms:
                                        X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
                                        mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2

                                        if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                                            X_to_pi_center_distance = Plevin.distance(X_pos_array, pi_center_array)

                                            if X_to_pi_center_distance < 6:
                                                XPCN_angle = Plevin.XPCN_angle(X_pos_array, pi_center_array, normal_vector)

                                                if XPCN_angle < 25.0:
                                                    alt_H = X_atom.altloc
                                                    H_atoms = neighbor_search.find_atoms(X_cra.atom.pos, alt=alt_H, min_dist=0.0, radius=1.3)
                                                    for H_atom in H_atoms:
                                                        H_cra, H_chain_name, H_residue_num, H_residue_name, H_element_name, H_atom_name, H_b_iso, H_pos, H_pos_array = get_info.atominfo(H_atom, model)
                                                        if H_cra.atom.element == gemmi.Element('H'):
                                                            XH_picenterAngle = Plevin.XH_picenterAngle(pi_center_array, X_pos_array, H_pos_array) 
                                                            
                                                            if XH_picenterAngle is not None and XH_picenterAngle > 120.0:                        
                                                                pi_residue_name = 'TRP_A'

                                                                Plevin_append_result.append_result(result, pdb_name, resolution, chain, residue, mean_b_factor,
                                                                    pi_residue_name, pi_center_array, normal_vector, pi_b_factor_mean, 
                                                                    X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_pos_array, X_b_iso,
                                                                    H_atom_name, H_pos_array, X_to_pi_center_distance, XH_picenterAngle, XPCN_angle)

        except Exception as e:
            print(f"Error processing file {filepath}: {e}")

        finally:
            pbar.update(1)

print(f"Detect {len(result)} XH-π interactions totally")

# save df to csv file
df = pd.DataFrame(result)
df.to_csv(os.path.join(output_dir, f'result_{method}.csv'), index=False)
print(f'The result has been saved to {os.path.join(output_dir, f"result_{method}.csv")}')
