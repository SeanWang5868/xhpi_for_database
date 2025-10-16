import numpy as np

def distance(pos1, pos2):
    distance = np.linalg.norm(pos1 - pos2)
    return distance

# Convert a string representation of a list to a numpy array
def string_to_array(s):
    return np.array([float(x) for x in s.strip('[]').split()])

# Calculate the angle between two vectors formed by three points: X, H, and center
def XH_picenterAngle(center, X, H):
    vector_HX = H - X
    vector_HC = H - center

    dot_product = np.dot(vector_HX, vector_HC)
    norm_HX = np.linalg.norm(vector_HX)
    norm_HC = np.linalg.norm(vector_HC)

    cos_theta = dot_product / (norm_HX * norm_HC)
    angle = np.arccos(cos_theta)  
    return np.degrees(angle)

# Calculate the angle between the vector X-pi_C and the normal vector
def XPCN_angle(X_pos, pi_center, pi_normal_vector):
    X = np.array(X_pos)
    Pi_C = np.array(pi_center)
    V = np.array(pi_normal_vector)
    
    # Calculate the vector XPi_C
    XPi_C = Pi_C - X
    
    # Calculate the dot product
    dot_product = np.dot(XPi_C, V)
    
    # Calculate the magnitudes
    norm_XPi_C = np.linalg.norm(XPi_C)
    norm_V = np.linalg.norm(V)
    
    # Calculate the cosine of the angle
    cos_theta = dot_product / (norm_XPi_C * norm_V)
    
    # Ensure the cosine value is within the valid range [-1, 1]
    cos_theta = np.clip(cos_theta, -1, 1)
    
    # Calculate the angle in radians and convert to degrees
    angle = np.arccos(cos_theta)
    angle_degrees = np.degrees(angle)
    
    # Adjust the angle to be <= 90 degrees if necessary
    if angle_degrees > 90:
        angle_degrees = 180 - angle_degrees
    
    return angle_degrees