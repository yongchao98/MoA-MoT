import math

def calculate_magnetic_field():
    """
    Calculates the longitudinal component of the magnetic field (B_z)
    at a specified target point due to a current element on a cylinder's surface.
    """
    # --- Step 1: Define Constants and Given Variables ---
    mu0 = 4 * math.pi * 1e-7  # Permeability of free space in H/m
    J_theta_mag = 40.0         # Magnitude of current density (interpreted as source strength)
    
    # Source point in cylindrical coordinates
    r_s = 0.34                 # meters
    theta_s_deg = 60.0         # degrees
    z_s = 0.5                  # meters
    
    # Target field point in Cartesian coordinates
    x_f = 0.1                  # meters
    y_f = -0.1                 # meters
    z_f = 0.2                  # meters

    # --- Step 2: Convert source point to Cartesian coordinates ---
    theta_s_rad = math.radians(theta_s_deg)
    x_s = r_s * math.cos(theta_s_rad)
    y_s = r_s * math.sin(theta_s_rad)
    
    # --- Step 3: Calculate necessary vectors in Cartesian coordinates ---
    # Displacement vector R = r_f - r_s
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    
    # Current density vector J (in the azimuthal direction a_theta)
    # a_theta = -sin(theta) * i + cos(theta) * j
    J_x = -J_theta_mag * math.sin(theta_s_rad)
    J_y = J_theta_mag * math.cos(theta_s_rad)
    J_z = 0.0

    # --- Step 4: Apply the Biot-Savart Law for the z-component ---
    # Numerator of the Biot-Savart law: z-component of (J x R)
    cross_product_z = J_x * R_y - J_y * R_x
    
    # Denominator of the Biot-Savart law: |R|^3
    R_mag_sq = R_x**2 + R_y**2 + R_z**2
    R_mag = math.sqrt(R_mag_sq)
    R_mag_cubed = R_mag**3
    
    # Prefactor
    mu0_over_4pi = mu0 / (4 * math.pi)
    
    # Final B_z calculation
    B_z = mu0_over_4pi * cross_product_z / R_mag_cubed

    # --- Step 5: Output the results, including the equation with numbers ---
    print("Calculation of the longitudinal component of the magnetic field (B_z).")
    
    print("\n--- Key Components for the Biot-Savart Law ---")
    print(f"Displacement Vector R = (x_f-x_s, y_f-y_s, z_f-z_s) = ({R_x:.4f}, {R_y:.4f}, {R_z:.4f}) m")
    print(f"Current Density Vector J = (-J_θ*sin(θ), J_θ*cos(θ), 0) = ({J_x:.4f}, {J_y:.4f}, {J_z:.4f}) A/m²")

    print("\n--- Final Equation: B_z = (μ₀/4π) * (J_x*R_y - J_y*R_x) / |R|³ ---")
    
    print("\nSubstituting the numerical values:")
    numerator_str = f"({J_x:.4f}) * ({R_y:.4f}) - ({J_y:.4f}) * ({R_x:.4f})"
    print(f"Numerator = J_x*R_y - J_y*R_x = {numerator_str} = {cross_product_z:.4f}")
    
    denominator_str = f"sqrt(({R_x:.4f})² + ({R_y:.4f})² + ({R_z:.4f})²)^3"
    print(f"Denominator = |R|³ = {denominator_str} = {R_mag_cubed:.4f}")
    
    print(f"\nB_z = ({mu0_over_4pi:.1e}) * ({cross_product_z:.4f}) / ({R_mag_cubed:.4f})")
    
    print("\n--- Final Answer ---")
    print(f"The magnitude of the longitudinal component of the magnetic field B_z is: {B_z:.4e} T")
    
    # Final raw value for the answer tag
    global final_answer
    final_answer = B_z
    
if __name__ == '__main__':
    final_answer = 0
    calculate_magnetic_field()
    print(f"<<<{final_answer}>>>")