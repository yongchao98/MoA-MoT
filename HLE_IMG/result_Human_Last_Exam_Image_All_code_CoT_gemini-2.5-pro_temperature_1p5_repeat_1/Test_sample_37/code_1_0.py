import math

def calculate_bz_density():
    """
    Calculates the longitudinal component of the magnetic field density
    at a target point due to a surface current element on a cylinder.
    """
    # --- 1. Define Constants and Given Values ---
    mu0_over_4pi = 1e-7  # H/m or T*m/A

    # Source point in cylindrical coordinates
    r_s_cyl = 0.34  # m
    theta_s_deg = 60  # degrees
    z_s = 0.5  # m

    # Target field point in Cartesian coordinates
    x_f = 0.1  # m
    y_f = -0.1  # m
    z_f = 0.2  # m
    
    # Magnitude of the azimuthal surface current density (K)
    # Assuming the given "J" in A/m^2 is a typo for K in A/m
    K_mag = 40.0  # A/m

    # --- 2. Convert Source Point to Cartesian Coordinates ---
    theta_s_rad = math.radians(theta_s_deg)
    x_s = r_s_cyl * math.cos(theta_s_rad)
    y_s = r_s_cyl * math.sin(theta_s_rad)
    # z_s is already given

    # --- 3. Calculate the Displacement Vector R ---
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    
    # --- 4. Calculate the Magnitude of R ---
    R_mag_sq = R_x**2 + R_y**2 + R_z**2
    R_mag = math.sqrt(R_mag_sq)
    R_mag_cubed = R_mag**3

    # --- 5. Determine the Current Density Vector K in Cartesian Coordinates ---
    # K is azimuthal, so its direction is the cylindrical theta_hat vector.
    # theta_hat = (-sin(theta), cos(theta), 0) in Cartesian coordinates.
    K_x = -K_mag * math.sin(theta_s_rad)
    K_y = K_mag * math.cos(theta_s_rad)
    # K_z = 0 for an azimuthal current

    # --- 6. Calculate the Z-Component of the Cross Product (K x R) ---
    K_cross_R_z = K_x * R_y - K_y * R_x
    
    # --- 7. Calculate the Magnetic Field Density Bz ---
    # b_z = (mu0/4pi) * (K x R)_z / R^3
    bz_density = (mu0_over_4pi * K_cross_R_z) / R_mag_cubed

    # --- 8. Print the calculation step-by-step ---
    print("--- Calculation of B_z Density ---")
    print(f"The formula used is: B_z_density = (μ₀/4π) * (K x R)_z / R³")
    print("\nIntermediate values:")
    print(f"  Source point (Cartesian): ({x_s:.4f}, {y_s:.4f}, {z_s:.4f}) m")
    print(f"  Target point (Cartesian): ({x_f:.4f}, {y_f:.4f}, {z_f:.4f}) m")
    print(f"  Displacement vector R: ({R_x:.4f}, {R_y:.4f}, {R_z:.4f}) m")
    print(f"  Current density vector K: ({K_x:.4f}, {K_y:.4f}, 0.0000) A/m")
    
    print("\nPlugging numbers into the formula:")
    print(f"  μ₀/4π = {mu0_over_4pi} T*m/A")
    print(f"  (K x R)_z = (K_x * R_y) - (K_y * R_x) = ({K_x:.4f} * {R_y:.4f}) - ({K_y:.4f} * {R_x:.4f}) = {K_cross_R_z:.4f} A*m")
    print(f"  R = sqrt({R_x:.4f}² + {R_y:.4f}² + {R_z:.4f}²) = {R_mag:.4f} m")
    print(f"  R³ = {R_mag:.4f}³ = {R_mag_cubed:.4f} m³")
    
    print(f"\nFinal Equation:")
    numerator = mu0_over_4pi * K_cross_R_z
    print(f"  B_z_density = ( {mu0_over_4pi} * {K_cross_R_z:.4f} ) / {R_mag_cubed:.4f}")
    print(f"  B_z_density = {numerator:.4e} / {R_mag_cubed:.4f}")
    
    print(f"\nResult:")
    print(f"  The magnitude of the longitudinal component of the magnetic field density B_z is: {bz_density:.4e} T/m²")
    
    return bz_density

# Execute the calculation and store the result
final_bz_density = calculate_bz_density()

# The final answer requested by the user prompt
# This is extracted from the printed output for the final return value
# Note: due to the ambiguity, the result is a field density in T/m^2.
# We present the numerical value.
# <<<1.2018e-05>>>