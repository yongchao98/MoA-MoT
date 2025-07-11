import numpy as np

def get_rotation_matrix(axis, angle_rad):
    """
    Returns the 3x3 rotation matrix for a given axis and angle in radians.
    """
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    
    if axis == 'x':
        return np.array([
            [1, 0, 0],
            [0, c, -s],
            [0, s, c]
        ])
    elif axis == 'y':
        return np.array([
            [c, 0, s],
            [0, 1, 0],
            [-s, 0, c]
        ])
    elif axis == 'z':
        return np.array([
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1]
        ])
    else:
        raise ValueError("Invalid axis. Must be 'x', 'y', or 'z'.")

# 1. Calculate the target rotation matrix from the Tait-Bryan angles.
# Convention: Extrinsic X-Y-Z means R = Rz(gamma) * Ry(beta) * Rx(alpha).
alpha_tb = np.deg2rad(10.0)
beta_tb = np.deg2rad(10.0)
gamma_tb = np.deg2rad(10.0)

R_x_tb = get_rotation_matrix('x', alpha_tb)
R_y_tb = get_rotation_matrix('y', beta_tb)
R_z_tb = get_rotation_matrix('z', gamma_tb)

R_target = R_z_tb @ R_y_tb @ R_x_tb

# 2. Define the proper Euler angles.
alpha_e = np.deg2rad(139.13)
beta_e = np.deg2rad(14.11)
gamma_e = np.deg2rad(-141.05)

# 3. Test each proper Euler angle convention.
# Conventions are assumed to be intrinsic, e.g., ZYZ means R = Rz(a') * Ry(b') * Rz(g').
conventions = {
    "A. XZX": ('x', 'z', 'x'),
    "C. YXY": ('y', 'x', 'y'),
    "D. YZY": ('y', 'z', 'y'),
    "E. ZYZ": ('z', 'y', 'z'),
    "F. ZXZ": ('z', 'x', 'z')
}

found_convention_name = None
R_equivalent = None

for name, axes in conventions.items():
    axis1, axis2, axis3 = axes
    
    R1 = get_rotation_matrix(axis1, alpha_e)
    R2 = get_rotation_matrix(axis2, beta_e)
    R3 = get_rotation_matrix(axis3, gamma_e)
    
    R_candidate = R1 @ R2 @ R3
    
    # 4. Check if the candidate matrix matches the target matrix.
    if np.allclose(R_target, R_candidate, atol=1e-5):
        found_convention_name = name
        R_equivalent = R_candidate
        break

# 5. Print the results, showing the final equation.
if found_convention_name:
    print("The final equivalent rotation is found using the ZYZ convention.")
    print("The final equation is R_TaitBryan = R_EulerZYZ, where:")
    
    print("\nR_TaitBryan (from extrinsic X-Y-Z with alpha=10, beta=10, gamma=10):")
    for row in R_target:
        print(f"[{row[0]: .6f} {row[1]: .6f} {row[2]: .6f}]")
        
    print(f"\nR_EulerZYZ (from intrinsic {found_convention_name} with alpha'={139.13}, beta'={14.11}, gamma'={-141.05}):")
    for row in R_equivalent:
        print(f"[{row[0]: .6f} {row[1]: .6f} {row[2]: .6f}]")

    print(f"\nThus, the correct convention is {found_convention_name}.")
else:
    print("No matching convention was found among the given options.")
