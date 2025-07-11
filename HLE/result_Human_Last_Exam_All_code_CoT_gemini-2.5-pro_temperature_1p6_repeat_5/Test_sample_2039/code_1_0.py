import numpy as np

def solve_rotation_convention():
    """
    This script determines the correct proper Euler angle convention that is equivalent
    to a given extrinsic Tait-Bryan rotation.
    """
    
    # --- Helper functions for creating rotation matrices from angles in degrees ---
    def R_x(theta_deg):
        """Rotation matrix around X-axis."""
        theta_rad = np.deg2rad(theta_deg)
        c, s = np.cos(theta_rad), np.sin(theta_rad)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    def R_y(theta_deg):
        """Rotation matrix around Y-axis."""
        theta_rad = np.deg2rad(theta_deg)
        c, s = np.cos(theta_rad), np.sin(theta_rad)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

    def R_z(theta_deg):
        """Rotation matrix around Z-axis."""
        theta_rad = np.deg2rad(theta_deg)
        c, s = np.cos(theta_rad), np.sin(theta_rad)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    # --- Step 1: Calculate the target rotation matrix from the Tait-Bryan angles ---
    alpha, beta, gamma = 10.0, 10.0, 10.0
    print(f"Calculating target rotation matrix for extrinsic X-Y-Z rotation with angles:")
    print(f"α = {alpha}°, β = {beta}°, γ = {gamma}°")
    # For an extrinsic X-Y-Z rotation, the matrix order is Rz * Ry * Rx
    print(f"Matrix equation: R_target = R_z({gamma}) @ R_y({beta}) @ R_x({alpha})")
    
    R_target = R_z(gamma) @ R_y(beta) @ R_x(alpha)
    print("\nTarget Rotation Matrix:")
    print(R_target)
    print("-" * 50)

    # --- Step 2: Define the Euler angles and conventions to test ---
    alpha_p, beta_p, gamma_p = 139.13, 14.11, -141.05
    print(f"Testing proper Euler angle conventions with angles:")
    print(f"α' = {alpha_p}°, β' = {beta_p}°, γ' = {gamma_p}°\n")

    conventions = {
        "A. XZX": lambda a, b, c: R_x(a) @ R_z(b) @ R_x(c),
        "B. XYZ": lambda a, b, c: R_x(a) @ R_y(b) @ R_z(c),
        "C. YXY": lambda a, b, c: R_y(a) @ R_x(b) @ R_y(c),
        "D. YZY": lambda a, b, c: R_y(a) @ R_z(b) @ R_y(c),
        "E. ZYZ": lambda a, b, c: R_z(a) @ R_y(b) @ R_z(c),
        "F. ZXZ": lambda a, b, c: R_z(a) @ R_x(b) @ R_z(c),
    }

    # --- Step 3 & 4: Iterate through conventions, compare matrices, and identify the match ---
    found_answer = False
    for name, func in conventions.items():
        if "XYZ" in name: continue # Skip XYZ as it's not a proper Euler convention

        R_test = func(alpha_p, beta_p, gamma_p)
        # Use a tolerance because the input Euler angles are rounded
        is_match = np.allclose(R_target, R_test, atol=1e-4)
        
        print(f"Testing Convention: {name}")
        if is_match:
            print("Result: MATCH FOUND")
            print("\nThe equivalent Euler angle representation is:")
            
            final_equation = ""
            if "XZX" in name:
                final_equation = f"Convention: XZX, Equation: R = R_x({alpha_p}) @ R_z({beta_p}) @ R_x({gamma_p})"
            elif "YXY" in name:
                final_equation = f"Convention: YXY, Equation: R = R_y({alpha_p}) @ R_x({beta_p}) @ R_y({gamma_p})"
            elif "YZY" in name:
                final_equation = f"Convention: YZY, Equation: R = R_y({alpha_p}) @ R_z({beta_p}) @ R_y({gamma_p})"
            elif "ZYZ" in name:
                final_equation = f"Convention: ZYZ, Equation: R = R_z({alpha_p}) @ R_y({beta_p}) @ R_z({gamma_p})"
            elif "ZXZ" in name:
                final_equation = f"Convention: ZXZ, Equation: R = R_z({alpha_p}) @ R_x({beta_p}) @ R_z({gamma_p})"
            
            print(final_equation)
            print("\nCalculated Matrix for this convention:")
            print(R_test)
            found_answer = True
            break
        else:
            print("Result: No Match")
        print("-" * 25)

    if not found_answer:
        print("No matching convention was found in the given options.")

solve_rotation_convention()