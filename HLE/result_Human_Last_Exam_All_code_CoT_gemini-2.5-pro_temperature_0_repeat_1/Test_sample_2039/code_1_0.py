import numpy as np

def get_rotation_matrix(axis, angle_deg):
    """Creates a rotation matrix for a given axis and angle in degrees."""
    theta = np.deg2rad(angle_deg)
    c, s = np.cos(theta), np.sin(theta)
    if axis.upper() == 'X':
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    elif axis.upper() == 'Y':
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    elif axis.upper() == 'Z':
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    else:
        raise ValueError("Axis must be 'X', 'Y', or 'Z'")

def main():
    """
    Finds the Euler angle convention that matches a given Tait-Bryan rotation.
    """
    # 1. Calculate the target rotation matrix from the extrinsic Tait-Bryan angles
    alpha, beta, gamma = 10.0, 10.0, 10.0
    
    Rx_alpha = get_rotation_matrix('X', alpha)
    Ry_beta = get_rotation_matrix('Y', beta)
    Rz_gamma = get_rotation_matrix('Z', gamma)
    
    # For extrinsic X-Y-Z, the order of multiplication is Rz * Ry * Rx
    R_target = Rz_gamma @ Ry_beta @ Rx_alpha
    
    print(f"The target rotation is given by an extrinsic rotation using X-Y-Z convention with angles:")
    print(f"alpha = {alpha} degrees")
    print(f"beta = {beta} degrees")
    print(f"gamma = {gamma} degrees\n")

    # 2. Define the candidate Euler angles and conventions
    alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05
    
    conventions = {
        "A": "XZX", "B": "XYZ", "C": "YXY",
        "D": "YZY", "E": "ZYZ", "F": "ZXZ"
    }
    
    print(f"We are checking for an equivalent intrinsic rotation with angles:")
    print(f"alpha' = {alpha_p} degrees")
    print(f"beta' = {beta_p} degrees")
    print(f"gamma' = {gamma_p} degrees\n")

    # 3. Iterate through conventions, calculate matrices, and compare
    correct_convention = None
    for key, conv_str in conventions.items():
        axis1, axis2, axis3 = conv_str[0], conv_str[1], conv_str[2]
        
        # For intrinsic rotations, the order of multiplication is R1 * R2 * R3
        R1 = get_rotation_matrix(axis1, alpha_p)
        R2 = get_rotation_matrix(axis2, beta_p)
        R3 = get_rotation_matrix(axis3, gamma_p)
        R_candidate = R1 @ R2 @ R3
        
        # Compare the candidate matrix with the target matrix
        if np.allclose(R_target, R_candidate, atol=1e-5):
            correct_convention = key
            print(f"Checking convention {key} ({conv_str}): Match found!")
            break
        else:
            print(f"Checking convention {key} ({conv_str}): No match.")

    if correct_convention:
        print(f"\nThe equivalent rotation is achieved using the {conventions[correct_convention]} convention.")
    else:
        print("\nNo matching convention was found among the options.")

if __name__ == "__main__":
    main()