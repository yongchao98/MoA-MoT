import numpy as np

def get_rotation_matrices():
    """Returns a dictionary of functions for elementary rotation matrices."""
    
    def Rx(theta_rad):
        c, s = np.cos(theta_rad), np.sin(theta_rad)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    def Ry(theta_rad):
        c, s = np.cos(theta_rad), np.sin(theta_rad)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

    def Rz(theta_rad):
        c, s = np.cos(theta_rad), np.sin(theta_rad)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        
    return {'X': Rx, 'Y': Ry, 'Z': Rz}

def main():
    """
    Finds the equivalent Euler angle convention for a given Tait-Bryan rotation.
    """
    # Tait-Bryan angles in degrees
    alpha_tb = 10.0
    beta_tb = 10.0
    gamma_tb = 10.0

    # Euler angles in degrees
    alpha_e = 139.13
    beta_e = 14.11
    gamma_e = -141.05

    # Convert all angles to radians
    angles_tb_rad = np.deg2rad([alpha_tb, beta_tb, gamma_tb])
    angles_e_rad = np.deg2rad([alpha_e, beta_e, gamma_e])
    
    # Get elementary rotation matrix functions
    R = get_rotation_matrices()

    # 1. Calculate the target rotation matrix for extrinsic X-Y-Z
    # The order of matrix multiplication is Rz * Ry * Rx
    target_matrix = R['Z'](angles_tb_rad[2]) @ R['Y'](angles_tb_rad[1]) @ R['X'](angles_tb_rad[0])

    # 2. Define candidate Euler conventions and test each one
    conventions = {
        'A': 'XZX',
        'C': 'YXY',
        'D': 'YZY',
        'E': 'ZYZ',
        'F': 'ZXZ'
    }

    correct_convention_key = None
    correct_convention_name = None

    for key, name in conventions.items():
        # Assuming intrinsic rotations, the multiplication order is left-to-right
        # e.g., for ZYZ: Rz(alpha') * Ry(beta') * Rz(gamma')
        axis1, axis2, axis3 = name[0], name[1], name[2]
        
        candidate_matrix = R[axis1](angles_e_rad[0]) @ R[axis2](angles_e_rad[1]) @ R[axis3](angles_e_rad[2])
        
        # 3. Compare matrices
        if np.allclose(target_matrix, candidate_matrix, atol=1e-4):
            correct_convention_key = key
            correct_convention_name = name
            break
            
    if correct_convention_name:
        print(f"The equivalent rotation is achieved using the '{correct_convention_name}' convention.")
        print("The final equivalent equation is:")
        print(f"Extrinsic_XYZ(alpha={alpha_tb}°, beta={beta_tb}°, gamma={gamma_tb}°) == "
              f"Intrinsic_{correct_convention_name}(alpha'={alpha_e}°, beta'={beta_e}°, gamma'={gamma_e}°)")
        print(f"\nFinal Answer Code: {correct_convention_key}")
    else:
        print("No matching convention was found among the choices.")

if __name__ == "__main__":
    main()