import numpy as np

def calculate_magnetic_field_component():
    """
    Calculates the longitudinal component of the magnetic field (B_z) at a target point
    due to a current element on the surface of a cylinder.
    """

    # --- Given constants and parameters ---
    mu_0 = 4 * np.pi * 1e-7  # Permeability of free space, H/m

    # Source point in cylindrical coordinates
    r_s_val = 0.34  # m
    theta_s_deg = 60  # degrees
    z_s_val = 0.5   # m

    # Target field point in Cartesian coordinates
    x_f_val = 0.1   # m
    y_f_val = -0.1  # m
    z_f_val = 0.2   # m

    # Azimuthal surface current density (assuming A/m, correcting the A/m^2 typo)
    K_theta_val = 40.0  # A/m

    # --- Step 1: Convert source point to Cartesian coordinates ---
    theta_s_rad = np.deg2rad(theta_s_deg)
    x_s_val = r_s_val * np.cos(theta_s_rad)
    y_s_val = r_s_val * np.sin(theta_s_rad)
    
    r_s_vec = np.array([x_s_val, y_s_val, z_s_val])
    r_f_vec = np.array([x_f_val, y_f_val, z_f_val])

    # --- Step 2: Calculate the displacement vector R and its magnitude ---
    R_vec = r_f_vec - r_s_vec
    R_mag = np.linalg.norm(R_vec)

    # --- Step 3: Define the surface current density vector K in Cartesian coordinates ---
    # The direction is azimuthal (theta_hat), which in Cartesian is (-sin(theta), cos(theta), 0)
    Kx_val = K_theta_val * (-np.sin(theta_s_rad))
    Ky_val = K_theta_val * (np.cos(theta_s_rad))
    Kz_val = 0.0
    K_vec = np.array([Kx_val, Ky_val, Kz_val])

    # --- Step 4: Calculate the cross product K x R ---
    cross_product_vec = np.cross(K_vec, R_vec)
    
    # --- Step 5: Calculate the B_z component ---
    # The factor (mu_0 / 4pi) is 1e-7
    mu0_over_4pi = 1e-7
    Bz_val = mu0_over_4pi * cross_product_vec[2] / (R_mag**3)

    # --- Print the detailed calculation steps ---
    print("--- Detailed Calculation ---")
    print(f"1. Source point (r', θ', z') = ({r_s_val} m, {theta_s_deg}°, {z_s_val} m)")
    print(f"   In Cartesian coordinates (x', y', z'):")
    print(f"   x' = {r_s_val} * cos({theta_s_deg}°) = {x_s_val:.5f} m")
    print(f"   y' = {r_s_val} * sin({theta_s_deg}°) = {y_s_val:.5f} m")
    print(f"   z' = {z_s_val:.5f} m")
    print("\n" + "-"*30 + "\n")
    
    print(f"2. Target point (x_f, y_f, z_f) = ({x_f_val} m, {y_f_val} m, {z_f_val} m)")
    print("\n" + "-"*30 + "\n")
    
    print(f"3. Displacement vector R = r_f - r' = (x_f-x', y_f-y', z_f-z')")
    print(f"   R_x = {x_f_val} - {x_s_val:.5f} = {R_vec[0]:.5f} m")
    print(f"   R_y = {y_f_val} - {y_s_val:.5f} = {R_vec[1]:.5f} m")
    print(f"   R_z = {z_f_val} - {z_s_val:.5f} = {R_vec[2]:.5f} m")
    print(f"   Magnitude |R| = sqrt({R_vec[0]:.5f}² + {R_vec[1]:.5f}² + {R_vec[2]:.5f}²) = {R_mag:.5f} m")
    print("\n" + "-"*30 + "\n")
    
    print(f"4. Surface current density vector K (at source point)")
    print(f"   K_θ = {K_theta_val} A/m")
    print(f"   K_x = {K_theta_val} * (-sin({theta_s_deg}°)) = {K_vec[0]:.5f} A/m")
    print(f"   K_y = {K_theta_val} * cos({theta_s_deg}°) = {K_vec[1]:.5f} A/m")
    print(f"   K_z = {K_vec[2]:.5f} A/m")
    print("\n" + "-"*30 + "\n")
    
    print(f"5. Cross product (K x R)_z = K_x*R_y - K_y*R_x")
    print(f"   (K x R)_z = ({K_vec[0]:.5f}) * ({R_vec[1]:.5f}) - ({K_vec[1]:.5f}) * ({R_vec[0]:.5f})")
    print(f"   (K x R)_z = {cross_product_vec[2]:.5f}")
    print("\n" + "-"*30 + "\n")
    
    print("6. Final Calculation for the longitudinal component B_z")
    print("   B_z = (μ₀ / 4π) * (K x R)_z / |R|³")
    print(f"   B_z = (1e-7) * ({cross_product_vec[2]:.5f}) / ({R_mag:.5f})³")
    print(f"   B_z = (1e-7) * ({cross_product_vec[2]:.5f}) / ({R_mag**3:.5f})")
    print(f"   B_z = {Bz_val:.4e} T")
    
    print("\n" + "="*30 + "\n")
    print(f"The magnitude of the longitudinal component B_z is {Bz_val:.4e} Tesla.")
    
    return Bz_val

if __name__ == '__main__':
    result = calculate_magnetic_field_component()
    # The final answer is wrapped in <<<>>> as requested.
    # print(f"\n<<<{result:.4e}>>>") # This would be for the final submission format.