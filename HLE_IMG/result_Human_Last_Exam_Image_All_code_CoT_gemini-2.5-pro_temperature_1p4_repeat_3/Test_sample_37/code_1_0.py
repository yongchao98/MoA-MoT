import math

def calculate_bz():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element at a source point on a cylinder's surface.
    """
    # --- Given Constants and Values ---
    mu_0 = 4 * math.pi * 1e-7  # Permeability of free space (H/m)
    # Assuming the given current density is a surface current density K_theta
    # with a typo in units.
    K_theta = 40.0  # Azimuthal surface current density (A/m)
    
    # Source point in cylindrical coordinates
    r_s = 0.34  # m
    theta_s_deg = 60.0  # degrees
    z_s = 0.5  # m
    
    # Target field point in Cartesian coordinates
    x_f = 0.1   # m
    y_f = -0.1  # m
    z_f = 0.2   # m
    
    # --- Step 1: Convert source coordinates to Cartesian ---
    theta_s_rad = math.radians(theta_s_deg)
    x_s = r_s * math.cos(theta_s_rad)
    y_s = r_s * math.sin(theta_s_rad)
    # z_s remains the same
    
    # --- Step 2: Calculate the displacement vector R and its magnitude ---
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    
    R_magnitude_sq = R_x**2 + R_y**2 + R_z**2
    R_magnitude = math.sqrt(R_magnitude_sq)
    R_magnitude_cubed = R_magnitude**3

    # --- Step 3: Express the current density vector K in Cartesian coordinates ---
    # K is in the azimuthal (theta_hat) direction.
    # theta_hat = -sin(theta) * i_hat + cos(theta) * j_hat
    K_x = -K_theta * math.sin(theta_s_rad)
    K_y = K_theta * math.cos(theta_s_rad)
    # K_z = 0
    
    # --- Step 4: Calculate the z-component of the cross product (K x R) ---
    cross_product_z = K_x * R_y - K_y * R_x
    
    # --- Step 5: Calculate Bz using the Biot-Savart Law ---
    # We are calculating dBz/dS, which gives the contribution per unit area.
    # To get a value in Tesla, we assume a source area dS = 1 m^2.
    # B_z = (mu_0 / (4 * pi)) * (K x R)_z / R^3
    mu_0_over_4pi = 1e-7
    B_z = mu_0_over_4pi * cross_product_z / R_magnitude_cubed

    # --- Print the breakdown of the calculation ---
    print("--- Calculation Breakdown ---")
    print(f"Source Point (Cartesian): (xs, ys, zs) = ({x_s:.4f}, {y_s:.4f}, {z_s:.4f}) m")
    print(f"Target Point (Cartesian): (xf, yf, zf) = ({x_f:.4f}, {y_f:.4f}, {z_f:.4f}) m")
    print("\nDisplacement Vector R = (xf-xs, yf-ys, zf-zs):")
    print(f"Rx = {x_f:.4f} - {x_s:.4f} = {R_x:.4f} m")
    print(f"Ry = {y_f:.4f} - {y_s:.4f} = {R_y:.4f} m")
    print(f"Rz = {z_f:.4f} - {z_s:.4f} = {R_z:.4f} m")
    print(f"Magnitude R = sqrt(Rx^2 + Ry^2 + Rz^2) = {R_magnitude:.4f} m")
    
    print("\nCurrent Density Vector K = (Kx, Ky, Kz):")
    print(f"Kx = -{K_theta} * sin({theta_s_deg}) = {K_x:.4f} A/m")
    print(f"Ky =  {K_theta} * cos({theta_s_deg}) = {K_y:.4f} A/m")
    print(f"Kz = 0.0 A/m")

    print("\nFinal Equation: Bz = (μ₀/4π) * (Kx*Ry - Ky*Rx) / R³")
    print(f"μ₀/4π = {mu_0_over_4pi:.2e} T·m/A")
    print(f"Numerator (Kx*Ry - Ky*Rx) = ({K_x:.4f} * {R_y:.4f}) - ({K_y:.4f} * {R_x:.4f}) = {cross_product_z:.4f} A")
    print(f"Denominator R³ = ({R_magnitude:.4f})³ = {R_magnitude_cubed:.4f} m³")
    
    print("\n--- Final Result ---")
    print(f"Bz = {mu_0_over_4pi:.2e} * {cross_product_z:.4f} / {R_magnitude_cubed:.4f}")
    print(f"The magnitude of the longitudinal component of the magnetic field Bz is: {abs(B_z):.4e} T")

if __name__ == '__main__':
    calculate_bz()