import numpy as np

def get_rotation_matrix(axis, angle_deg):
    """
    Computes the 3x3 rotation matrix for a given axis and angle.
    :param axis: 'X', 'Y', or 'Z'
    :param angle_deg: The rotation angle in degrees.
    :return: 3x3 numpy array representing the rotation matrix.
    """
    theta_rad = np.deg2rad(angle_deg)
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    
    if axis.upper() == 'X':
        return np.array([[1, 0, 0],
                         [0, c, -s],
                         [0, s, c]])
    elif axis.upper() == 'Y':
        return np.array([[c, 0, s],
                         [0, 1, 0],
                         [-s, 0, c]])
    elif axis.upper() == 'Z':
        return np.array([[c, -s, 0],
                         [s, c, 0],
                         [0, 0, 1]])
    else:
        raise ValueError("Axis must be 'X', 'Y', or 'Z'")

# --- Step 1: Calculate the reference rotation matrix ---

# Initial Tait-Bryan angles in degrees
alpha, beta, gamma = 10, 10, 10

# For an extrinsic X-Y-Z rotation, the matrix is Rz * Ry * Rx
R_ref = get_rotation_matrix('Z', gamma) @ get_rotation_matrix('Y', beta) @ get_rotation_matrix('X', alpha)

print("The reference rotation matrix from the X-Y-Z convention with angles ({}, {}, {}) is:".format(alpha, beta, gamma))
print(R_ref)
print("\n" + "="*70 + "\n")

# --- Step 2: Test each Euler angle convention ---

# Equivalent proper Euler angles in degrees
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05

# Define the conventions to be tested
conventions = {
    'A': 'XZX',
    'B': 'XYZ',
    'C': 'YXY',
    'D': 'YZY',
    'E': 'ZYZ',
    'F': 'ZXZ'
}

# Iterate through the conventions and find the correct one
correct_choice = None
correct_convention_name = None

for choice, name in conventions.items():
    # Construct the rotation matrix for the current convention
    R1 = get_rotation_matrix(name[0], alpha_p)
    R2 = get_rotation_matrix(name[1], beta_p)
    R3 = get_rotation_matrix(name[2], gamma_p)
    R_test = R1 @ R2 @ R3
    
    # Check if the test matrix is close to the reference matrix
    if np.allclose(R_ref, R_test, atol=1e-4):
        correct_choice = choice
        correct_convention_name = name
        break

# --- Step 3: Print the final result and verification ---

print(f"The equivalent Euler angles are alpha'={alpha_p}, beta'={beta_p}, gamma'={gamma_p}")

if correct_choice:
    print(f"The matching convention is '{correct_convention_name}', which is choice {correct_choice}.")

    # Re-calculate the correct matrix to show the equation
    R1 = get_rotation_matrix(correct_convention_name[0], alpha_p)
    R2 = get_rotation_matrix(correct_convention_name[1], beta_p)
    R3 = get_rotation_matrix(correct_convention_name[2], gamma_p)
    R_test_correct = R1 @ R2 @ R3

    print(f"\nVerification:")
    print(f"The rotation matrix for the {correct_convention_name} convention with these angles is:")
    print(R_test_correct)
    
    print("\nAs we can see, the two matrices are virtually identical, confirming the result.")
else:
    print("No matching convention was found among the choices.")

<<<E>>>