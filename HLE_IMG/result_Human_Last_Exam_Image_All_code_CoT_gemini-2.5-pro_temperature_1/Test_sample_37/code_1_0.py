import math

def calculate_bz():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element on a cylinder's surface.

    Assumptions made to solve the problem:
    1. The given "azimuthal current density (J)" of 40 A/m^2 is interpreted as
       a surface current density K with a value of 40 A/m. The unit A/m^2 is
       assumed to be a typo.
    2. To calculate a finite magnetic field value (B) in Tesla from a single source
       point, we calculate the contribution from a patch of area 1 m^2 (ΔS=1).
    """

    # --- Given constants and parameters ---
    mu0_over_4pi = 1e-7  # H/m
    
    # Source point in cylindrical coordinates
    r_s = 0.34  # m
    theta_s_deg = 60  # degrees
    z_s = 0.5   # m
    
    # Target field point in Cartesian coordinates
    x_f = 0.1   # m
    y_f = -0.1  # m
    z_f = 0.2   # m
    
    # Azimuthal surface current density (K_theta)
    K_theta = 40  # A/m (based on assumption 1)
    
    # Area element
    delta_S = 1.0 # m^2 (based on assumption 2)

    # --- Step-by-step calculation ---

    # 1. Convert source point coordinates to Cartesian
    theta_s_rad = math.radians(theta_s_deg)
    x_s = r_s * math.cos(theta_s_rad)
    y_s = r_s * math.sin(theta_s_rad)
    # z_s is already in Cartesian

    # 2. Calculate the displacement vector R = r_f - r_s
    Rx = x_f - x_s
    Ry = y_f - y_s
    Rz = z_f - z_s

    # 3. Calculate the magnitude of the displacement vector R
    R_squared = Rx**2 + Ry**2 + Rz**2
    R = math.sqrt(R_squared)

    # 4. Express the azimuthal current density vector K in Cartesian coordinates
    # K = K_theta * a_theta, where a_theta = [-sin(theta), cos(theta), 0]
    Kx = -K_theta * math.sin(theta_s_rad)
    Ky = K_theta * math.cos(theta_s_rad)
    Kz = 0

    # 5. Calculate the z-component of the cross product (K x R)
    # (K x R)_z = Kx * Ry - Ky * Rx
    Cz = Kx * Ry - Ky * Rx

    # 6. Calculate the longitudinal component of the magnetic field (Bz)
    Bz = mu0_over_4pi * Cz / (R**3) * delta_S
    
    # --- Print the results and the final equation ---
    print("This calculation is based on the Biot-Savart law for a surface current element:")
    print("B_z = (μ₀ / 4π) * [ (K × R)_z / |R|³ ] * ΔS\n")
    
    print("--- Input Values ---")
    print(f"Target Point (x_f, y_f, z_f) = ({x_f}, {y_f}, {z_f}) m")
    print(f"Source Point (r_s, θ_s, z_s) = ({r_s} m, {theta_s_deg}°, {z_s} m)")
    print(f"Surface Current Density K_θ = {K_theta} A/m")
    print(f"Permeability constant μ₀/4π = {mu0_over_4pi} H/m\n")

    print("--- Intermediate Calculations ---")
    print(f"Displacement Vector R = (x_f - r_s*cos(θ_s), y_f - r_s*sin(θ_s), z_f - z_s)")
    print(f"R = ({x_f} - {x_s:.4f}, {y_f} - {y_s:.4f}, {z_f} - {z_s})")
    print(f"R = ({Rx:.4f}, {Ry:.4f}, {Rz:.4f}) m")
    print(f"Magnitude |R| = sqrt({Rx:.4f}² + {Ry:.4f}² + {Rz:.4f}²) = {R:.4f} m")
    print(f"Z-component of Cross Product (K × R)_z = K_x*R_y - K_y*R_x")
    print(f"(K × R)_z = ({Kx:.4f})*({Ry:.4f}) - ({Ky:.4f})*({Rx:.4f}) = {Cz:.4f} A·m\n")

    print("--- Final Equation with Numbers ---")
    # Using a simplified expression for Cz: Cz = K_theta * (r_s - x_f*cos(theta) - y_f*sin(theta))
    # This avoids calculating Kx, Ky explicitly in the printout for clarity.
    Cz_simplified_part = r_s - x_f * math.cos(theta_s_rad) - y_f * math.sin(theta_s_rad)
    
    print(f"B_z = ({mu0_over_4pi}) * ( {K_theta} * ( {r_s} - {x_f}*cos({theta_s_deg}°) - ({y_f})*sin({theta_s_deg}°) ) ) / ({R:.4f})³")
    print(f"B_z = ({mu0_over_4pi}) * ( {K_theta} * ({Cz_simplified_part:.4f}) ) / ({R**3:.4f})")
    print(f"B_z = ({mu0_over_4pi}) * ( {Cz:.4f} ) / ({R**3:.4f})")
    print(f"B_z = {Bz:.4e} T")
    
    return Bz

# Execute the function and capture the final answer
final_answer = calculate_bz()
print(f"\n<<<The magnitude of the longitudinal component of the magnetic field (Bz) is approximately {final_answer:.4e} T.>>>")
print(f"<<<{final_answer}>>>")
