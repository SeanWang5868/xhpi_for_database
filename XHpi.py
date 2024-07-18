import os
import sys
import config
import pandas as pd
import gzip
import gemmi
from tqdm import tqdm
from datetime import datetime
import wget
from calculate import get_info, calculator, append_result, rename, printer, Plevin, Plevin_append_result
from addH import addH

printer.print_xhpi()

# Create a folder with today's date as the name
today = datetime.now().strftime('%Y-%m-%d')
output_dir = os.path.join(".", today)
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

# Get user input for the method
method = input("Select the detection system (Hudson or Plevin): ").strip().lower()
if method not in ['hudson', 'plevin']:
    print("Invalid method. Please enter 'Hudson' or 'Plevin'.")
    sys.exit(1)

# Process the downloaded files and calculate XHPI
gz_files = [os.path.join(output_dir, f"{pdb_name}.cif.gz") for pdb_name in pdb_names]

addH.main()

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
            pdb_name = os.path.basename(filepath).replace('.cif.gz', '')  # Get the PDB name

            # Read the gzipped CIF file
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
                                # Find X neighbor atoms
                                alt_pi = pi_atoms[0].altloc
                                X_atoms = neighbor_search.find_atoms(pi_center, alt=alt_pi, min_dist=0.0, radius=6)

                                if method == 'hudson':
                                    # Calculate using Hudson method
                                    for X_atom in X_atoms:
                                        X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
                                        mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2

                                        if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                                            X_to_pi_center_distance = calculator.distance(X_pos_array, pi_center_array)

                                            if X_to_pi_center_distance <= 4.5:
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
                                    # Calculate using Plevin method
                                    for X_atom in X_atoms:
                                        X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
                                        mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2

                                        if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                                            X_to_pi_center_distance = Plevin.distance(X_pos_array, pi_center_array)

                                            if X_to_pi_center_distance < 4.3:
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
                                # Find X neighbor atoms
                                alt_pi = pi_atoms[0].altloc
                                X_atoms = neighbor_search.find_atoms(pi_center, alt=alt_pi, min_dist=0.0, radius=6)

                                if method == 'hudson':
                                    # Calculate using Hudson method
                                    for X_atom in X_atoms:
                                        X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
                                        mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2

                                        if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                                            X_to_pi_center_distance = calculator.distance(X_pos_array, pi_center_array)

                                            if X_to_pi_center_distance <= 4.5:
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
                                    # Calculate using Plevin method
                                    for X_atom in X_atoms:
                                        X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
                                        mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2

                                        if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                                            X_to_pi_center_distance = Plevin.distance(X_pos_array, pi_center_array)

                                            if X_to_pi_center_distance < 4.3:
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

# Save the calculation results to a CSV file
df = pd.DataFrame(result)
df.to_csv(os.path.join(output_dir, f'result_{method}.csv'), index=False)
print(f'The result has been saved to {os.path.join(output_dir, f"result_{method}.csv")}')

