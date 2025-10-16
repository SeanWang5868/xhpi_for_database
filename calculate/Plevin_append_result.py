
def append_result(result, pdb_name, resolution, chain, residue, mean_b_factor,
                  pi_residue_name, pi_center_array, normal_vector, pi_b_factor_mean, 
                  X_chain_name, X_residue_num, X_residue_name, X_element_name, X_atom_name, X_pos_array, X_b_iso,
                  H_atom_name, H_pos_array, 
                  X_to_pi_center_distance, XH_picenterAngle, XPCN_angle):

# def append_result(result, pdb_name, resolution, pi_residue_name, residue, chain, cra, X_b_factor, target_atom_pos, hydrogen_atom_name, hydrogen_atom_pos, 
#                   X_to_pi_center_distance, angle_in_degrees, projection_distance, X_residue_chain_name):
    
    result.append({
        'pdb_name': pdb_name,
        'resolution': resolution,
        'mean_b_factor': mean_b_factor,

        'pi_chain_name': chain.name,
        'pi_residue_number': residue.seqid.num,
        'pi_residue_name': pi_residue_name,
        'pi_center_array': pi_center_array,
        'pi_normal_vector': normal_vector,
        'pi_b_factor': pi_b_factor_mean,

        'X_chain_name': X_chain_name,
        # 'X_residue_id':X_residue_id,
        'X_residue_num': X_residue_num,
        'X_residue_name': X_residue_name,
        'X_element_name': X_element_name,
        'X_atom_name': X_atom_name,
        'X_pos': X_pos_array,
        'X_b_factor': X_b_iso,
        'H_atom_name': H_atom_name,
        'H_pos': H_pos_array,

        'X_to_pi_center_distance': X_to_pi_center_distance,
        'XH_picenterAngle': XH_picenterAngle,
        'XPCN_angle': XPCN_angle,
        # 'interaction_energy': interaction_energy
    })