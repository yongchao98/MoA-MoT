import numpy as np

def solve_rotation_convention():
    """
    Calculates the target rotation matrix from Tait-Bryan angles and
    finds the matching proper Euler angle convention from a list of candidates.
    """
    # -- Step 1: Calculate the target rotation matrix from Tait-Bryan angles --

    # Tait-Bryan angles in degrees
    alpha_tb_deg, beta_tb_deg, gamma_tb_deg = 10.0, 10.0, 10.0
    
    # Convert angles to radians for trigonometric functions
    alpha_tb = np.deg2rad(alpha_tb_deg)
    beta_tb = np.deg2rad(beta_tb_deg)
    gamma_tb = np.deg2rad(gamma_tb_deg)

    # Define the basic rotation matrices
    def Rx(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    def Ry(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

    def Rz(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    # Calculate the target rotation matrix for extrinsic X-Y-Z convention
    # R_target = R_Z(gamma) * R_Y(beta) * R_X(alpha)
    R_target = Rz(gamma_tb) @ Ry(beta_tb) @ Rx(alpha_tb)

    print("--- Task ---")
    print(f"Find the Euler angle convention for angles α'={139.13}°, β'={14.11}°, γ'={-141.05}°")
    print(f"that is equivalent to a Tait-Bryan X-Y-Z rotation with angles α={alpha_tb_deg}°, β={beta_tb_deg}°, γ={gamma_tb_deg}°.")
    print("\n--- Calculation ---")
    print("1. Target Rotation Matrix (from Tait-Bryan angles):")
    print(np.round(R_target, 4))
    
    # -- Step 2: Test each candidate Euler angle convention --

    # Given equivalent proper Euler angles in degrees
    alpha_e_deg, beta_e_deg, gamma_e_deg = 139.13, 14.11, -141.05

    # Convert to radians
    alpha_e = np.deg2rad(alpha_e_deg)
    beta_e = np.deg2rad(beta_e_deg)
    gamma_e = np.deg2rad(gamma_e_deg)

    # Define the conventions to test
    # Key: Answer Choice, Value: (Name, Rotation Function Sequence)
    conventions = {
        'A': ('XZX', (Rx, Rz, Rx)),
        'B': ('XYZ', (Rx, Ry, Rz)),
        'C': ('YXY', (Ry, Rx, Ry)),
        'D': ('YZY', (Ry, Rz, Ry)),
        'E': ('ZYZ', (Rz, Ry, Rz)),
        'F': ('ZXZ', (Rz, Rx, Rz)),
    }

    print("\n2. Testing Candidate Conventions:")
    found_match = False
    result_key = None
    
    for key, (name, (R1_func, R2_func, R3_func)) in conventions.items():
        # Calculate the candidate rotation matrix for the current convention
        R_candidate = R1_func(alpha_e) @ R2_func(beta_e) @ R3_func(gamma_e)
        
        # Compare with the target matrix using a tolerance for floating point errors
        # and rounding in the problem's given angles.
        is_match = np.allclose(R_target, R_candidate, atol=1e-3)
        
        print(f"\n- Convention {name} (Choice {key}):")
        print(np.round(R_candidate, 4))
        if is_match:
            print("  --> Match found!")
            found_match = True
            result_key = key
        else:
            print("  --> No match.")
    
    print("\n--- Conclusion ---")
    if found_match:
        print(f"The correct convention is {conventions[result_key][0]}.")
        # Returning the answer key in the required format
        print(f"<<<{result_key}>>>")
    else:
        print("No matching convention was found among the options.")

solve_rotation_convention()