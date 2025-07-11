import numpy as np

def get_rotation_matrix(angle_deg, axis):
    """Creates a rotation matrix for a given angle in degrees and axis."""
    angle_rad = np.deg2rad(angle_deg)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    if axis.upper() == 'X':
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    elif axis.upper() == 'Y':
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    elif axis.upper() == 'Z':
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    else:
        raise ValueError("Axis must be 'X', 'Y', or 'Z'")

# Step 1: Calculate the target rotation matrix from the Tait-Bryan angles
alpha_tb, beta_tb, gamma_tb = 10.0, 10.0, 10.0

# For extrinsic XYZ, the matrix is Rz(gamma) * Ry(beta) * Rx(alpha)
Rx_tb = get_rotation_matrix(alpha_tb, 'X')
Ry_tb = get_rotation_matrix(beta_tb, 'Y')
Rz_tb = get_rotation_matrix(gamma_tb, 'Z')
R_target = Rz_tb @ Ry_tb @ Rx_tb

print(f"Original Tait-Bryan rotation: X({alpha_tb}°) Y({beta_tb}°) Z({gamma_tb}°)")
print("Target Rotation Matrix:")
print(R_target)
print("-" * 30)

# Step 2: Define the Euler angles to test
alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05

# Step 3 & 4: Test each convention from the answer choices
conventions = {
    'A': 'XZX',
    'B': 'XYZ', # Tait-Bryan, not proper Euler, but we'll check
    'C': 'YXY',
    'D': 'YZY',
    'E': 'ZYZ',
    'F': 'ZXZ'
}

found_convention = None
final_R = None
final_R1, final_R2, final_R3 = None, None, None

for key, conv_str in conventions.items():
    axis1, axis2, axis3 = conv_str[0], conv_str[1], conv_str[2]
    
    # For intrinsic rotations, the matrix order matches the convention string
    # e.g., ZYZ is Rz(a) * Ry(b) * Rz(g)
    # For Tait-Bryan (like XYZ), it's Rz(g) * Ry(b) * Rx(a)
    if conv_str in ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX']: # Tait-Bryan
         R1 = get_rotation_matrix(alpha_p, axis1)
         R2 = get_rotation_matrix(beta_p, axis2)
         R3 = get_rotation_matrix(gamma_p, axis3)
         R_test = R3 @ R2 @ R1
    else: # Proper Euler
        R1 = get_rotation_matrix(alpha_p, axis1)
        R2 = get_rotation_matrix(beta_p, axis2)
        R3 = get_rotation_matrix(gamma_p, axis3)
        R_test = R1 @ R2 @ R3

    # Compare with the target matrix
    if np.allclose(R_target, R_test, atol=1e-4):
        found_convention = f"{key}. {conv_str}"
        final_R = R_test
        final_R1, final_R2, final_R3 = R1, R2, R3
        break

# Step 5: Output the final result
if found_convention:
    print(f"Found matching convention: {found_convention}")
    print(f"The equivalent rotation is given by the {found_convention.split('. ')[1]} convention with angles:")
    print(f"α' = {alpha_p}°, β' = {beta_p}°, γ' = {gamma_p}°\n")
    
    axis1, axis2, axis3 = found_convention.split('. ')[1]
    
    print("The final equation is:")
    if found_convention.split('. ')[1] in ['XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX']:
        print(f"R = R_{axis3}({gamma_p}°) @ R_{axis2}({beta_p}°) @ R_{axis1}({alpha_p}°)")
        print("\nEach matrix in the equation is:")
        print(f"R_{axis3}({gamma_p}°) =\n{final_R3}\n")
        print(f"R_{axis2}({beta_p}°) =\n{final_R2}\n")
        print(f"R_{axis1}({alpha_p}°) =\n{final_R1}\n")
    else:
        print(f"R = R_{axis1}({alpha_p}°) @ R_{axis2}({beta_p}°) @ R_{axis3}({gamma_p}°)")
        print("\nEach matrix in the equation is:")
        print(f"R_{axis1}({alpha_p}°) =\n{final_R1}\n")
        print(f"R_{axis2}({beta_p}°) =\n{final_R2}\n")
        print(f"R_{axis3}({gamma_p}°) =\n{final_R3}\n")

    print("Resulting Matrix R =")
    print(final_R)
    print("\nThis matches the target matrix.")
    
    answer_key = found_convention.split('.')[0]
    print(f"\n<<<E>>>")

else:
    print("No matching convention was found among the choices.")
