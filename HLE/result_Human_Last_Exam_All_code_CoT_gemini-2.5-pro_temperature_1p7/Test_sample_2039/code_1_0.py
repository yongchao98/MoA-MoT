import numpy as np

def get_rotation_matrices():
    """Returns functions to generate rotation matrices for X, Y, Z axes."""
    def Rx(theta_deg):
        theta_rad = np.deg2rad(theta_deg)
        c, s = np.cos(theta_rad), np.sin(theta_rad)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    def Ry(theta_deg):
        theta_rad = np.deg2rad(theta_deg)
        c, s = np.cos(theta_rad), np.sin(theta_rad)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

    def Rz(theta_deg):
        theta_rad = np.deg2rad(theta_deg)
        c, s = np.cos(theta_deg), np.sin(theta_deg)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return Rx, Ry, Rz

def main():
    """
    Finds the Euler angle convention that matches a given Tait-Bryan rotation.
    """
    # Get rotation matrix functions
    Rx, Ry, Rz = get_rotation_matrices()

    # --- Step 1: Calculate the Target Rotation Matrix ---
    alpha, beta, gamma = 10.0, 10.0, 10.0
    # For extrinsic X-Y-Z rotation, the matrix is Rz(gamma) * Ry(beta) * Rx(alpha)
    R_target = Rz(gamma) @ Ry(beta) @ Rx(alpha)

    print("The rotation can be described by the Tait-Bryan angles alpha=10, beta=10, gamma=10 with the extrinsic X-Y-Z convention.")
    print(f"The corresponding equation is R = Rz({gamma}) @ Ry({beta}) @ Rx({alpha})")
    print("The resulting target rotation matrix is:")
    print(R_target)
    print("-" * 30)

    # --- Step 2 & 3: Iterate through conventions and find the match ---
    alpha_prime, beta_prime, gamma_prime = 139.13, 14.11, -141.05

    conventions = {
        "A": {"name": "XZX", "funcs": [Rx, Rz, Rx]},
        "B": {"name": "XYZ", "funcs": [Rx, Ry, Rz]},
        "C": {"name": "YXY", "funcs": [Ry, Rx, Ry]},
        "D": {"name": "YZY", "funcs": [Ry, Rz, Ry]},
        "E": {"name": "ZYZ", "funcs": [Rz, Ry, Rz]},
        "F": {"name": "ZXZ", "funcs": [Rz, Rx, Rz]},
    }
    
    correct_choice = None
    
    for choice, props in conventions.items():
        name = props["name"]
        R1, R2, R3 = props["funcs"]
        
        # Calculate the candidate matrix for the current convention
        R_candidate = R1(alpha_prime) @ R2(beta_prime) @ R3(gamma_prime)
        
        # Check if the candidate matrix matches the target matrix
        if np.allclose(R_target, R_candidate, atol=1e-4):
            correct_choice = choice
            print(f"Found a match with the {name} convention.")
            print("The final equation using this convention is:")
            print(f"R = {name[0]}({alpha_prime}) @ {name[1]}({beta_prime}) @ {name[2]}({gamma_prime})")
            print("The resulting matrix is:")
            print(R_candidate)
            print("\nComparing the two matrices confirms they are equivalent.")
            break

    # This part of the code is for the final answer submission format.
    if correct_choice:
        # We need to explicitly print the formatted answer string
        # to ensure it's captured correctly.
        final_answer_string = f"<<<{correct_choice}>>>"
        print(f"\nFinal Answer: {final_answer_string}")


if __name__ == "__main__":
    main()