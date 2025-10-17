import gemmi
from calculate import get_info, rename, Plevin_calculator, Plevin_append_result

def detect_plevin(pdb_name, resolution, model, chain, structure, residue, atoms_dict):
    """
    Detects potential pi-interactions for a given residue.
    """
    result_temp = []

    # Prepare for neighbor search
    neighbor_search = gemmi.NeighborSearch(structure[0], structure.cell, 8)
    neighbor_search.populate(include_h=True)
    pi_atoms = [atom for atom in residue if atom.name in atoms_dict[residue.name]]

    # Ensure all atoms of the pi-system are present
    if len(pi_atoms) == len(atoms_dict[residue.name]):
        pi_center, pi_center_array, pi_atoms_array, normal_vector, pi_b_factor_mean = get_info.piinfo(pi_atoms)
        alt_pi = pi_atoms[0].altloc

        # Find neighboring atoms (potential X atoms)
        X_atoms = neighbor_search.find_atoms(pi_center, alt=alt_pi, min_dist=0.0, radius=8)

        # Evaluate each potential interacting atom
        for X_atom in X_atoms:
            X_cra, X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_b_iso, X_pos, X_pos_array = get_info.atominfo(X_atom, model)
            
            # Check if atom type is C, N, O, or S
            if X_cra.atom.element in [gemmi.Element('C'), gemmi.Element('N'), gemmi.Element('O'), gemmi.Element('S')]:
                X_to_pi_center_distance = Plevin_calculator.distance(X_pos_array, pi_center_array)

                # First distance cutoff
                if X_to_pi_center_distance < 8.0:
                    XPCN_angle = Plevin_calculator.XPCN_angle(X_pos_array, pi_center_array, normal_vector)

                    # First angle cutoff
                    if XPCN_angle < 45.0:
                        alt_H = X_atom.altloc
                        H_atoms = neighbor_search.find_atoms(X_cra.atom.pos, alt=alt_H, min_dist=0.0, radius=1.3)
                        
                        for H_atom in H_atoms:
                            H_cra, _, _, _, _, H_atom_name, _, _, H_pos_array = get_info.atominfo(H_atom, model)
                            
                            # Ensure the neighbor is a Hydrogen atom
                            if H_cra.atom.element == gemmi.Element('H'):
                                XH_picenterAngle = Plevin_calculator.XH_picenterAngle(pi_center_array, X_pos_array, H_pos_array)
                                
                                # Second angle cutoff
                                if XH_picenterAngle is not None and XH_picenterAngle > 90.0:
                                    
                                    # --- START: Added Logic for 'xhpi' Column ---
                                    # Check if the interaction meets the specific geometric criteria.
                                    if (X_to_pi_center_distance < 4.3 and
                                        XH_picenterAngle > 120.0 and
                                        XPCN_angle < 25.0):
                                        xhpi = 1
                                    else:
                                        xhpi = 0
                                    # ---- END: Added Logic for 'xhpi' Column ----

                                    mean_b_factor = (pi_b_factor_mean + X_b_iso) / 2
                                    pi_residue_name = rename.pi_rename(residue, len(atoms_dict[residue.name]))
                                    
                                    # Append the result, now including the xhpi_flag
                                    # NOTE: Ensure your Plevin_append_result.append_result function
                                    # is updated to accept this new 'xhpi_flag' argument.
                                    Plevin_append_result.append_result(
                                        result_temp, pdb_name, resolution, chain, residue, mean_b_factor,
                                        pi_residue_name, pi_center_array, normal_vector, pi_b_factor_mean, 
                                        X_chain_name, X_residue_num, X_residue_name, X_element_name, 
                                        X_atom_name, X_pos_array, X_b_iso, H_atom_name, H_pos_array, 
                                        X_to_pi_center_distance, XH_picenterAngle, XPCN_angle,
                                        xhpi  # Pass the new flag here
                                    )
    return result_temp