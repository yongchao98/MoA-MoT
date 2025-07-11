import math

def calculate_bz():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element on a cylinder's surface.
    """
    # --- Given Constants and Values ---
    mu0_over_4pi = 1e-7  # mu0 / (4*pi) in H/m
    
    # Target field point (Cartesian)
    xf, yf, zf = 0.1, -0.1, 0.2  # meters
    
    # Source point (Cylindrical)
    r_s = 0.34  # meters
    theta_s_deg = 60.0  # degrees
    z_s = 0.5  # meters
    
    # Azimuthal surface current density (assuming K_theta)
    K_theta = 40.0  # A/m
    
    # Assumed differential surface area element
    dS = 1.0 # m^2

    # --- Step 1: Convert source point to Cartesian coordinates ---
    theta_s_rad = math.radians(theta_s_deg)
    xs = r_s * math.cos(theta_s_rad)
    ys = r_s * math.sin(theta_s_rad)
    # zs is already given
    
    # --- Step 2: Calculate the displacement vector R and its magnitude ---
    Rx = xf - xs
    Ry = yf - ys
    Rz = zf - z_s
    R_vec = (Rx, Ry, Rz)
    
    R_mag_sq = Rx**2 + Ry**2 + Rz**2
    R_mag = math.sqrt(R_mag_sq)
    
    # --- Step 3: Express the current density vector K in Cartesian coordinates ---
    # The azimuthal unit vector in Cartesian is (-sin(theta), cos(theta), 0)
    Kx = -K_theta * math.sin(theta_s_rad)
    Ky = K_theta * math.cos(theta_s_rad)
    Kz = 0
    K_vec = (Kx, Ky, Kz)
    
    # --- Step 4: Calculate the z-component of the cross product (K x R) ---
    cross_product_z = Kx * Ry - Ky * Rx
    
    # --- Step 5: Calculate Bz using the Biot-Savart Law ---
    Bz = mu0_over_4pi * (cross_product_z / (R_mag**3)) * dS

    # --- Print the results and the final equation ---
    print("--- Problem Setup ---")
    print(f"Target Point (xf, yf, zf) = ({xf}, {yf}, {zf}) m")
    print(f"Source Point (r, theta, z) = ({r_s}, {theta_s_deg} deg, {z_s}) m")
    print(f"Azimuthal Surface Current Density K_theta = {K_theta} A/m")
    print(f"Permeability constant (mu_0 / 4pi) = {mu0_over_4pi} H/m\n")
    
    print("--- Calculation Steps ---")
    print(f"1. Source point in Cartesian (xs, ys, zs) = ({xs:.4f}, {ys:.4f}, {zs:.4f}) m")
    print(f"2. Displacement vector R = (xf-xs, yf-ys, zf-zs) = ({Rx:.4f}, {Ry:.4f}, {Rz:.4f}) m")
    print(f"3. Magnitude of R = {R_mag:.4f} m")
    print(f"4. Current density vector K = (Kx, Ky, Kz) = ({Kx:.4f}, {Ky:.4f}, {Kz:.4f}) A/m")
    print(f"5. Z-component of cross product (K x R)_z = {cross_product_z:.4f} A\n")

    print("--- Final Equation ---")
    print("B_z = (mu_0 / 4*pi) * (K x R)_z / R^3 * dS")
    print(f"B_z = ({mu0_over_4pi}) * ({cross_product_z:.4f}) / ({R_mag:.4f})^3 * {dS}")
    print(f"B_z = ({mu0_over_4pi}) * ({cross_product_z:.4f}) / ({R_mag**3:.4f}) * {dS}\n")

    print(f"--- Final Answer ---")
    print(f"The magnitude of the longitudinal component of the magnetic field B_z is: {Bz:.4e} T")
    
    # For automated checking
    # print(f"<<<{Bz}>>>")

calculate_bz()