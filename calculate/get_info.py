import numpy as np
import gemmi

def atominfo(atom, structure):
    cra = atom.to_cra(structure)

    chain_name = cra.chain.name,
    residue_id = cra.residue.entity_id,
    residue_num = cra.residue.seqid.num,
    residue_name = cra.residue.name,
    element_name = cra.atom.element.name,
    atom_name = cra.atom.name,
    b_iso = cra.atom.b_iso,
    pos = cra.atom.pos,
    pos_array = np.array([cra.atom.pos.x, cra.atom.pos.y, cra.atom.pos.z])
    
    return cra, chain_name[0], residue_num[0], residue_name[0], element_name[0], atom_name[0], b_iso[0], pos, pos_array

def piinfo(pi_atoms):
    coords = np.array([atom.pos for atom in pi_atoms])
    center = np.mean(coords, axis=0)
    pi_center = gemmi.Position(center[0], center[1], center[2])
    pi_center_array = np.array([pi_center.x, pi_center.y, pi_center.z])

    pi_atoms_array = np.array([[atom.pos.x, atom.pos.y, atom.pos.z] for atom in pi_atoms])
    v1 = pi_atoms_array[1] - pi_atoms_array[0]
    v2 = pi_atoms_array[2] - pi_atoms_array[0]
    normal_vector = np.cross(v1, v2)
    normal_vector = normal_vector / np.linalg.norm(normal_vector)  # 归一化法向量

    pi_b_factor = []
    for atom in pi_atoms: # b_factor
        pi_b_factor.append(atom.b_iso)
        pi_b_factor_mean = np.mean(pi_b_factor)

    return pi_center, pi_center_array, pi_atoms_array, normal_vector, pi_b_factor_mean