import numpy as np

def get_rotation_matrix(axis, angle_deg):
    """
    Returns the 3x3 rotation matrix for a given axis and angle.
    """
    theta = np.radians(angle_deg)
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
    Main function to find the equivalent Euler angle convention.
    """
    # 1. Calculate the reference rotation matrix from the extrinsic Tait-Bryan angles
    alpha_tb, beta_tb, gamma_tb = 10.0, 10.0, 10.0
    
    # For extrinsic X-Y-Z, the matrix is Rz(gamma) * Ry(beta) * Rx(alpha)
    R_x_tb = get_rotation_matrix('X', alpha_tb)
    R_y_tb = get_rotation_matrix('Y', beta_tb)
    R_z_tb = get_rotation_matrix('Z', gamma_tb)
    
    R_ref = R_z_tb @ R_y_tb @ R_x_tb

    # 2. Define the equivalent proper Euler angles
    alpha_e, beta_e, gamma_e = 139.13, 14.11, -141.05

    # 3. Define the candidate conventions (assuming intrinsic rotations)
    # For intrinsic A-B-C, the matrix is RA(alpha) * RB(beta) * RC(gamma)
    conventions = {
        "A. XZX": (('X', alpha_e), ('Z', beta_e), ('X', gamma_e)),
        "B. XYZ": (('X', alpha_e), ('Y', beta_e), ('Z', gamma_e)),
        "C. YXY": (('Y', alpha_e), ('X', beta_e), ('Y', gamma_e)),
        "D. YZY": (('Y', alpha_e), ('Z', beta_e), ('Y', gamma_e)),
        "E. ZYZ": (('Z', alpha_e), ('Y', beta_e), ('Z', gamma_e)),
        "F. ZXZ": (('Z', alpha_e), ('X', beta_e), ('Z', gamma_e))
    }
    
    correct_convention_label = None
    final_equation_angles = None

    # 4. Iterate through conventions, calculate matrix, and compare
    for label, rotations in conventions.items():
        axis1, angle1 = rotations[0]
        axis2, angle2 = rotations[1]
        axis3, angle3 = rotations[2]

        R1 = get_rotation_matrix(axis1, angle1)
        R2 = get_rotation_matrix(axis2, angle2)
        R3 = get_rotation_matrix(axis3, angle3)
        
        R_cand = R1 @ R2 @ R3
        
        # Check if the candidate matrix is close to the reference matrix
        if np.allclose(R_ref, R_cand, atol=1e-3):
            correct_convention_label = label
            final_equation_angles = (axis1, angle1, axis2, angle2, axis3, angle3)
            final_matrix = R_cand
            break

    # Print the final result including the equation
    if correct_convention_label:
        print(f"The initial rotation is given by extrinsic X-Y-Z with angles ({alpha_tb:.2f}°, {beta_tb:.2f}°, {gamma_tb:.2f}°).")
        print("This results in the reference matrix:")
        print(R_ref)
        print("\n" + "="*50 + "\n")
        
        print(f"The equivalent rotation is found using the '{correct_convention_label}' convention.")
        
        axis1, angle1, axis2, angle2, axis3, angle3 = final_equation_angles
        print("The final equation for the equivalent rotation matrix is:")
        print(f"R = R{axis1.lower()}({angle1:.2f}°) * R{axis2.lower()}({angle2:.2f}°) * R{axis3.lower()}({angle3:.2f}°)")
        
        print("\nSubstituting the values, each number in the equation is:")
        print(f"Angle 1: {angle1:.2f}°")
        print(f"Angle 2: {angle2:.2f}°")
        print(f"Angle 3: {angle3:.2f}°")
        
        print("\nThis calculation yields the matrix:")
        print(final_matrix)
        
        print("\nThis matrix matches the reference matrix, so the correct convention is", correct_convention_label[0])
    else:
        print("No matching convention was found among the options.")
        
    # Return the final answer in the required format
    if correct_convention_label:
        print(f"\n<<<{correct_convention_label[0]}>>>")

if __name__ == '__main__':
    main()