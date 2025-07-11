import numpy as np

def solve_rotation_convention():
    """
    Finds the proper Euler angle convention that matches a given Tait-Bryan rotation.
    """

    # --- Step 1: Calculate the Reference Rotation Matrix (Tait-Bryan) ---
    # Tait-Bryan angles in degrees for an extrinsic X-Y-Z rotation
    alpha_deg, beta_deg, gamma_deg = 10.0, 10.0, 10.0
    
    # Convert angles to radians for numpy trigonometric functions
    alpha = np.radians(alpha_deg)
    beta = np.radians(beta_deg)
    gamma = np.radians(gamma_deg)

    # Define the basic rotation matrices
    def Rx(theta):
        return np.array([
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)]
        ])

    def Ry(theta):
        return np.array([
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)]
        ])

    def Rz(theta):
        return np.array([
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1]
        ])

    # For extrinsic X-Y-Z rotation, the matrix is R = Rz(gamma) * Ry(beta) * Rx(alpha)
    R_ref = Rz(gamma) @ Ry(beta) @ Rx(alpha)
    print("Initial Rotation: Tait-Bryan X(10.0)Y(10.0)Z(10.0)")
    
    # --- Step 2 & 3: Test Each Proper Euler Angle Convention ---
    # Equivalent Euler angles in degrees
    alpha_p_deg, beta_p_deg, gamma_p_deg = 139.13, 14.11, -141.05
    
    # Convert angles to radians
    alpha_p = np.radians(alpha_p_deg)
    beta_p = np.radians(beta_p_deg)
    gamma_p = np.radians(gamma_p_deg)
    
    # Dictionary of conventions and their matrix multiplication order (intrinsic)
    conventions = {
        "A. XZX": Rx(alpha_p) @ Rz(beta_p) @ Rx(gamma_p),
        "C. YXY": Ry(alpha_p) @ Rx(beta_p) @ Ry(gamma_p),
        "D. YZY": Ry(alpha_p) @ Rz(beta_p) @ Ry(gamma_p),
        "E. ZYZ": Rz(alpha_p) @ Ry(beta_p) @ Rz(gamma_p),
        "F. ZXZ": Rz(alpha_p) @ Rx(beta_p) @ Rz(gamma_p),
    }

    # --- Step 4: Compare Matrices and Find the Match ---
    found_match = False
    for name, R_candidate in conventions.items():
        if np.allclose(R_ref, R_candidate, atol=1e-3):
            print(f"\nFound a match: Convention {name}")
            print(f"The rotation from this convention is equivalent to the initial rotation using angles:")
            print(f"alpha' = {alpha_p_deg} degrees")
            print(f"beta'  = {beta_p_deg} degrees")
            print(f"gamma' = {gamma_p_deg} degrees")
            
            # Print the final matrices to show they are equivalent
            np.set_printoptions(precision=4, suppress=True)
            print("\nReference Tait-Bryan Matrix:")
            print(R_ref)
            print(f"\nEquivalent {name} Matrix:")
            print(R_candidate)
            
            found_match = True
            break
            
    if not found_match:
        print("\nNo matching convention was found among the options.")

# Run the solver
solve_rotation_convention()