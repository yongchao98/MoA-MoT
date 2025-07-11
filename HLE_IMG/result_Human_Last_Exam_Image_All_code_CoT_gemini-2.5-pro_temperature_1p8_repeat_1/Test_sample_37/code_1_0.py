import math

def calculate_magnetic_field_component():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element at a source point on a cylinder's surface.
    """
    # --- 1. Define constants and given values ---
    # The constant mu_0 / (4*pi)
    mu0_over_4pi = 1e-7  # T*m/A or H/m

    # Azimuthal surface current density in A/m.
    # The problem statement's unit A/m^2 is assumed to be a typo for A/m.
    K_theta_magnitude = 40.0

    # Source point in cylindrical coordinates (r, theta, z)
    r_s_cyl = 0.34
    theta_s_deg = 60.0
    z_s = 0.5

    # Target field point in Cartesian coordinates (x, y, z)
    x_f = 0.1
    y_f = -0.1
    z_f = 0.2

    # --- 2. Convert source point to Cartesian coordinates ---
    theta_s_rad = math.radians(theta_s_deg)
    x_s = r_s_cyl * math.cos(theta_s_rad)
    y_s = r_s_cyl * math.sin(theta_s_rad)

    # --- 3. Calculate the displacement vector R = r_f - r_s and its magnitude ---
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    R_mag = math.sqrt(R_x**2 + R_y**2 + R_z**2)

    # --- 4. Define the surface current density vector K in Cartesian coordinates ---
    # The azimuthal unit vector theta_hat = -sin(theta) * i_hat + cos(theta) * j_hat
    K_x = -K_theta_magnitude * math.sin(theta_s_rad)
    K_y = K_theta_magnitude * math.cos(theta_s_rad)
    # K_z is 0 as the current is purely azimuthal

    # --- 5. Calculate the z-component of the cross product (K x R) ---
    cross_product_z = K_x * R_y - K_y * R_x

    # --- 6. Calculate the z-component of the magnetic field ---
    Bz = mu0_over_4pi * cross_product_z / (R_mag**3)

    # --- 7. Print the results showing the numbers in the final equation ---
    print("This script calculates the longitudinal magnetic field component B_z.")
    print("Based on the formula: B_z = (μ₀/4π) * (K x R)_z / R³\n")
    print("--- Key Values ---")
    print(f"Displacement Vector R = ({R_x:.4f}, {R_y:.4f}, {R_z:.4f}) m")
    print(f"Current Density Vector K = ({K_x:.4f}, {K_y:.4f}, 0.0000) A/m")
    
    print("\n--- Equation Components ---")
    print(f"Numerator (K x R)_z = ({K_x:.4f}) * ({R_y:.4f}) - ({K_y:.4f}) * ({R_x:.4f}) = {cross_product_z:.4f} A")
    print(f"Denominator R³ = ({R_mag:.4f} m)³ = {R_mag**3:.4f} m³")
    print(f"Constant (μ₀/4π) = {mu0_over_4pi:.2e} H/m")
    
    print("\n--- Final Equation ---")
    print(f"B_z = ({mu0_over_4pi:.2e}) * ({cross_product_z:.4f}) / ({R_mag**3:.4f})")
    print(f"B_z = {Bz:.4e} T")


if __name__ == "__main__":
    calculate_magnetic_field_component()