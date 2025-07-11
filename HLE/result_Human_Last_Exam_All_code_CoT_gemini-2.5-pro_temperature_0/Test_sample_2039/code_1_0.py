import numpy as np

def get_rotation_matrix(alpha_deg, beta_deg, gamma_deg, convention):
    """
    Calculates the rotation matrix for a given convention and angles.
    Angles are applied in the order alpha, beta, gamma.
    The resulting matrix is for an extrinsic rotation sequence.
    """
    # Convert angles from degrees to radians
    alpha = np.radians(alpha_deg)
    beta = np.radians(beta_deg)
    gamma = np.radians(gamma_deg)

    # Basic rotation matrices
    Rx_a = np.array([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
    Ry_a = np.array([[np.cos(alpha), 0, np.sin(alpha)], [0, 1, 0], [-np.sin(alpha), 0, np.cos(alpha)]])
    Rz_a = np.array([[np.cos(alpha), -np.sin(alpha), 0], [np.sin(alpha), np.cos(alpha), 0], [0, 0, 1]])

    Rx_b = np.array([[1, 0, 0], [0, np.cos(beta), -np.sin(beta)], [0, np.sin(beta), np.cos(beta)]])
    Ry_b = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
    Rz_b = np.array([[np.cos(beta), -np.sin(beta), 0], [np.sin(beta), np.cos(beta), 0], [0, 0, 1]])

    Rx_g = np.array([[1, 0, 0], [0, np.cos(gamma), -np.sin(gamma)], [0, np.sin(gamma), np.cos(gamma)]])
    Ry_g = np.array([[np.cos(gamma), 0, np.sin(gamma)], [0, 1, 0], [-np.sin(gamma), 0, np.cos(gamma)]])
    Rz_g = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])

    # Extrinsic rotation composition: R = R_gamma * R_beta * R_alpha
    if convention == 'XYZ':
        return Rz_g @ Ry_b @ Rx_a
    elif convention == 'XZX':
        return Rx_g @ Rz_b @ Rx_a
    elif convention == 'YXY':
        return Ry_g @ Rx_b @ Ry_a
    elif convention == 'YZY':
        return Ry_g @ Rz_b @ Ry_a
    elif convention == 'ZYZ':
        return Rz_g @ Ry_b @ Rz_a
    elif convention == 'ZXZ':
        return Rz_g @ Rx_b @ Rz_a
    else:
        raise ValueError("Unknown convention")

# --- 1. Calculate the reference rotation matrix ---
# Initial Tait-Bryan angles
alpha_tb = 10.0
beta_tb = 10.0
gamma_tb = 10.0
# The convention is extrinsic X-Y-Z
R_ref = get_rotation_matrix(alpha_tb, beta_tb, gamma_tb, 'XYZ')

# --- 2. Test each candidate convention ---
# Given proper Euler angles
alpha_p = 139.13
beta_p = 14.11
gamma_p = -141.05

# List of candidate conventions from the answer choices
# Note: XYZ is Tait-Bryan, not proper Euler, but we include it for completeness.
candidate_conventions = {
    "A": "XZX",
    "B": "XYZ",
    "C": "YXY",
    "D": "YZY",
    "E": "ZYZ",
    "F": "ZXZ"
}

# --- 3. Compare and find the match ---
found = False
for key, convention_name in candidate_conventions.items():
    # Calculate the rotation matrix for the candidate convention
    R_candidate = get_rotation_matrix(alpha_p, beta_p, gamma_p, convention_name)
    
    # Check if the candidate matrix is close to the reference matrix
    if np.allclose(R_ref, R_candidate, atol=1e-4):
        print(f"The correct convention is {key}. {convention_name}.")
        print(f"This convention uses the angles:")
        print(f"alpha' = {alpha_p} degrees")
        print(f"beta'  = {beta_p} degrees")
        print(f"gamma' = {gamma_p} degrees")
        print("This produces the same rotation as the initial Tait-Bryan rotation.")
        found = True
        break

if not found:
    print("No matching convention was found among the options.")
