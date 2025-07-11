import math

def calculate_magnetic_field_z():
    """
    Calculates the longitudinal component of the magnetic field (B_z) at a target point
    due to a current density element on a cylinder surface.
    """
    # --- Step 1: Define Given Parameters ---
    # Permeability of free space (H/m)
    mu_0 = 4 * math.pi * 1e-7

    # Target field point (Cartesian coordinates, m)
    x_f, y_f, z_f = 0.1, -0.1, 0.2

    # Source point (Cylindrical coordinates)
    r_s = 0.34      # m
    theta_s_deg = 60  # degrees
    z_s = 0.5       # m

    # Azimuthal current density (A/m^2)
    J_theta = 40.0

    print("--- Input Parameters ---")
    print(f"Permeability of free space (μ₀): {mu_0:.4e} H/m")
    print(f"Target field point (x_f, y_f, z_f): ({x_f}, {y_f}, {z_f}) m")
    print(f"Source point (r, θ, z): ({r_s} m, {theta_s_deg}°, {z_s} m)")
    print(f"Azimuthal current density (J_θ): {J_theta} A/m²\n")

    # --- Step 2: Coordinate Transformation ---
    print("--- Step 1: Coordinate Transformation ---")
    theta_s_rad = math.radians(theta_s_deg)
    
    # Convert source point to Cartesian coordinates
    x_s = r_s * math.cos(theta_s_rad)
    y_s = r_s * math.sin(theta_s_rad)
    print(f"Source point in Cartesian (x_s, y_s, z_s): ({x_s:.4f}, {y_s:.4f}, {z_s:.1f}) m")

    # Convert current density vector J to Cartesian coordinates
    # J = J_theta * a_theta = J_theta * (-sin(theta) i + cos(theta) j)
    J_x = J_theta * (-math.sin(theta_s_rad))
    J_y = J_theta * (math.cos(theta_s_rad))
    J_z = 0.0 # Azimuthal current has no z-component
    print(f"Current density vector in Cartesian (J_x, J_y, J_z): ({J_x:.4f}, {J_y:.4f}, {J_z:.1f}) A/m²\n")

    # --- Step 3: Calculate Displacement Vector R ---
    print("--- Step 2: Calculate Displacement Vector R ---")
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    print(f"Displacement vector R = (x_f-x_s, y_f-y_s, z_f-z_s)")
    print(f"R = ({R_x:.4f}, {R_y:.4f}, {R_z:.4f}) m\n")
    
    # --- Step 4: Apply Biot-Savart Law for B_z ---
    print("--- Step 3: Apply Biot-Savart Law ---")
    # Magnitude of R
    R_mag_sq = R_x**2 + R_y**2 + R_z**2
    R_mag = math.sqrt(R_mag_sq)
    R_mag_cubed = R_mag**3
    
    # z-component of the cross product (J x R)
    cross_product_z = J_x * R_y - J_y * R_x
    print(f"(J x R)_z = J_x*R_y - J_y*R_x = ({J_x:.4f})*({R_y:.4f}) - ({J_y:.4f})*({R_x:.4f}) = {cross_product_z:.4f}")

    # Calculate B_z assuming dS = 1 m^2
    # B_z = (μ₀ / 4π) * (J x R)_z / |R|^3
    B_z = (mu_0 / (4 * math.pi)) * cross_product_z / R_mag_cubed
    
    print("\nFinal Equation for B_z:")
    print("B_z = (μ₀ / 4π) * (J_x*R_y - J_y*R_x) / (R_x² + R_y² + R_z²)^(3/2)")
    print("Substituting the values:")
    print(f"B_z = ({mu_0 / (4 * math.pi):.1e}) * (({J_x:.4f})*({R_y:.4f}) - ({J_y:.4f})*({R_x:.4f})) / (({R_x:.4f})² + ({R_y:.4f})² + ({R_z:.4f})²)^(3/2)")
    print(f"B_z = ({mu_0 / (4 * math.pi):.1e}) * ({cross_product_z:.4f}) / ({R_mag_sq:.4f})^(1.5)")
    print(f"B_z = ({mu_0 / (4 * math.pi):.1e}) * ({cross_product_z:.4f}) / {R_mag_cubed:.4f}\n")
    
    print("--- Result ---")
    print(f"The magnitude of the longitudinal component of the magnetic field B_z is: {B_z:.4e} T")
    
    return B_z

if __name__ == '__main__':
    final_bz = calculate_magnetic_field_z()
    # The final answer in a simple format as requested
    # Expressing in microtesla (uT) for readability, 1 T = 1e6 uT
    final_bz_uT = final_bz * 1e6
    # print(f"\nFinal Answer: {final_bz_uT:.3f} uT")
    # print(f"<<<{final_bz_uT:.3f}>>>") # As microtesla
    print(f"<<<{final_bz:.4e}>>>") # As Tesla in scientific notation
