import math

def calculate_bz():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element on a cylinder surface based on the Biot-Savart law.
    """
    # --- 1. Given Information ---
    # Permeability of free space
    mu_0 = 4 * math.pi * 1e-7  # H/m
    
    # Target field point in Cartesian coordinates (m)
    xf, yf, zf = 0.1, -0.1, 0.2
    r_f = (xf, yf, zf)
    
    # Source point in Cylindrical coordinates (m, degrees, m)
    r_s_cyl, theta_s_deg, z_s_cyl = 0.34, 60, 0.5
    
    # Azimuthal surface current density magnitude (A/m^2, interpreted as A/m for calculation)
    J_theta_mag = 40.0

    print("--- Step-by-step Calculation ---")

    # --- 2. Coordinate Conversion ---
    # Convert angle from degrees to radians
    theta_s_rad = math.radians(theta_s_deg)
    
    # Convert source point to Cartesian coordinates
    xs = r_s_cyl * math.cos(theta_s_rad)
    ys = r_s_cyl * math.sin(theta_s_rad)
    zs = z_s_cyl
    r_s = (xs, ys, zs)
    
    print(f"\n1. Source Point (r_s):")
    print(f"   Cylindrical: (r={r_s_cyl}, theta={theta_s_deg}°, z={z_s_cyl}) m")
    print(f"   Cartesian:   (x_s={xs:.4f}, y_s={ys:.4f}, z_s={zs:.4f}) m")

    print(f"\n2. Target Point (r_f):")
    print(f"   Cartesian:   (x_f={xf}, y_f={yf}, z_f={zf}) m")

    # --- 3. Determine Key Vectors ---
    # Displacement vector R = r_f - r_s
    Rx = xf - xs
    Ry = yf - ys
    Rz = zf - zs
    R = (Rx, Ry, Rz)
    
    # Magnitude of R
    R_mag = math.sqrt(Rx**2 + Ry**2 + Rz**2)
    
    print(f"\n3. Displacement Vector (R = r_f - r_s):")
    print(f"   R = ({Rx:.4f}, {Ry:.4f}, {Rz:.4f}) m")
    print(f"   Magnitude |R| = {R_mag:.4f} m")

    # Current density vector J in Cartesian coordinates
    # J is in the azimuthal direction a_theta = (-sin(theta), cos(theta), 0)
    Jx = J_theta_mag * -math.sin(theta_s_rad)
    Jy = J_theta_mag * math.cos(theta_s_rad)
    Jz = 0
    J = (Jx, Jy, Jz)
    
    print(f"\n4. Current Density Vector (J):")
    print(f"   J (Cartesian) = ({Jx:.4f}, {Jy:.4f}, {Jz:.4f}) A/m")

    # --- 4. Apply Biot-Savart Law ---
    # We need the z-component of the cross product C = J x R
    # Cz = Jx * Ry - Jy * Rx
    Cz = Jx * Ry - Jy * Rx
    
    print(f"\n5. Z-component of Cross Product (C_z = J_x*R_y - J_y*R_x):")
    print(f"   C_z = ({Jx:.4f}) * ({Ry:.4f}) - ({Jy:.4f}) * ({Rx:.4f})")
    print(f"   C_z = {Cz:.4f}")

    # Calculate Bz using the Biot-Savart law formula
    # Bz = (mu_0 / 4*pi) * Cz / |R|^3
    constant = mu_0 / (4 * math.pi) # This is 1e-7
    Bz = constant * Cz / (R_mag**3)
    
    print(f"\n6. Final Calculation for B_z:")
    print(f"   B_z = (μ₀ / 4π) * C_z / |R|³")
    print(f"   B_z = ({constant:.1e}) * ({Cz:.4f}) / ({R_mag:.4f})³")
    print(f"   B_z = ({constant:.1e}) * ({Cz:.4f}) / ({R_mag**3:.4f})")
    print(f"\n--- Result ---")
    print(f"The magnitude of the longitudinal component of the magnetic field B_z is: {abs(Bz):.4e} T")
    
    return abs(Bz)

# Run the calculation and print the final answer in the required format
final_answer = calculate_bz()
# The final answer is formatted as requested.
# For example: <<<9.8>>> or <<<C>>>
# The calculated value is approximately 1.2015e-5
# Let's format it as requested.
print(f"\n<<<{final_answer:.4e}>>>")
