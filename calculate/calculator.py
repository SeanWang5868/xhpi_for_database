import gemmi
import numpy as np

import gemmi  # Assuming gemmi package is correctly imported

def distance(pos1, pos2):
    distance = np.linalg.norm(pos1 - pos2)
    return distance

def theta_angle(pi_center_array, X_pos_array, H_pos_array, normal_vector):
    
    # the vector from X to the pi center and the XH vector
    X_pi_vector = pi_center_array - X_pos_array
    XH_vector = H_pos_array - X_pos_array

    # Calculate the vector projection length
    projection_length = np.dot(XH_vector, X_pi_vector) / np.linalg.norm(X_pi_vector)

    # XPi_D = np.linalg.norm(X_pos_array - pi_center_array)
    # HPi_D = np.linalg.norm(H_pos_array - pi_center_array)
    
    # Determine whether XH bond points to the pi ring
    if projection_length > 0:
    # if XPi_D > HPi_D:
        cosine_angle = np.dot(normal_vector, XH_vector) / (np.linalg.norm(normal_vector) * np.linalg.norm(XH_vector))
        angle_in_degrees = np.degrees(np.arccos(cosine_angle))
        if angle_in_degrees >= 90:
            angle_in_degrees = 180 - angle_in_degrees
        return angle_in_degrees
    return None

# def projection_distance(normal_vector, pi_atoms, pi_center, c_atom):
#     pi_center = np.array(pi_center)
#     c_atom = np.array(c_atom)

#     d = np.dot(pi_atoms[0], normal_vector)
#     t = (d - np.dot(c_atom, normal_vector)) / np.dot(normal_vector, normal_vector)

#     projection = c_atom + t * normal_vector
#     distance = np.linalg.norm(projection - pi_center)
#     return distance

def projection_distance(normal_vector, pi_center, c_atom):
    # pi_center = np.array(pi_center)
    # c_atom = np.array(c_atom)

    t = np.dot(normal_vector, pi_center - c_atom) / np.dot(normal_vector, normal_vector)

    projection = c_atom + t * normal_vector
    distance = np.linalg.norm(projection - pi_center)
    return distance


def interaction_strength(X_atom_pos, H_atom_pos, X_to_pi_center_distance, th_angle, projection_distance):
    # Calculate CH bond length
    CH_bond_length = np.linalg.norm(H_atom_pos - X_atom_pos)
    
    # Normalize distances if necessary (example scaling)
    X_to_pi_center_distance_norm = X_to_pi_center_distance / 10.0
    projection_distance_norm = projection_distance / 5.0
    
    # Calculate interaction strength using a weighted combination
    interaction_strength = (1 / X_to_pi_center_distance_norm * 0.5 + 
                            1 / th_angle * 0.2 + 
                            1 / projection_distance_norm * 0.2 +
                            1 / CH_bond_length * 0.1)
    
    return interaction_strength

