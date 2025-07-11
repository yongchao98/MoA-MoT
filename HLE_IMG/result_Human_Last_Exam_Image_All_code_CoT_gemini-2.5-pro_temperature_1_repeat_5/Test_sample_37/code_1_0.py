import math

def calculate_magnetic_field_component():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current density at a source point on a cylinder.
    """
    # --- GIVEN VALUES ---
    # Permeability of free space constant
    mu0_over_4pi = 1e-7  # H/m

    # Target field point in Cartesian coordinates (m)
    xf, yf, zf = 0.1, -0.1, 0.2

    # Source point in Cylindrical coordinates (m, degrees, m)
    rs, theta_deg, zs = 0.34, 60.0, 0.5

    # Magnitude of the azimuthal current density (A/m^2)
    J_theta = 40.0

    # --- STEP 1: Convert source point to Cartesian coordinates ---
    theta_rad = math.radians(theta_deg)
    xs = rs * math.cos(theta_rad)
    ys = rs * math.sin(theta_rad)

    # --- STEP 2: Calculate the displacement vector R ---
    Rx = xf - xs
    Ry = yf - ys
    Rz = zf - zs

    # --- STEP 3: Convert the current density vector J to Cartesian coordinates ---
    # The azimuthal unit vector theta_hat in Cartesian is (-sin(theta), cos(theta), 0)
    Jx = -J_theta * math.sin(theta_rad)
    Jy = J_theta * math.cos(theta_rad)
    # Jz is 0 for an azimuthal current

    # --- STEP 4: Calculate the components of the Bz formula ---
    # Numerator: z-component of the cross product (J x R)
    numerator = Jx * Ry - Jy * Rx
    
    # Denominator: Magnitude of R cubed
    R_mag_squared = Rx**2 + Ry**2 + Rz**2
    R_mag = math.sqrt(R_mag_squared)
    denominator = R_mag**3
    
    # Final calculation for Bz
    Bz = mu0_over_4pi * numerator / denominator
    
    # --- OUTPUT RESULTS ---
    print("This calculation determines the longitudinal magnetic field component (Bz) using the Biot-Savart law.")
    print("The formula is: Bz = (mu0 / 4pi) * (Jx * Ry - Jy * Rx) / (Rx^2 + Ry^2 + Rz^2)^(3/2)\n")

    print("--- Values plugged into the equation ---")
    print(f"mu0 / 4pi = {mu0_over_4pi} H/m")
    print(f"Jx = {Jx:.5f} A/m^2")
    print(f"Jy = {Jy:.5f} A/m^2")
    print(f"Rx = {Rx:.5f} m")
    print(f"Ry = {Ry:.5f} m")
    print(f"Rz = {Rz:.5f} m\n")

    print("--- Step-by-step evaluation of the equation ---")
    print(f"Bz = {mu0_over_4pi} * (({Jx:.4f}) * ({Ry:.4f}) - ({Jy:.4f}) * ({Rx:.4f})) / (({Rx:.4f})^2 + ({Ry:.4f})^2 + ({Rz:.4f})^2)^(1.5)")
    print(f"Bz = {mu0_over_4pi} * ({numerator:.4f}) / ({R_mag_squared:.4f})^(1.5)")
    print(f"Bz = {mu0_over_4pi} * ({numerator:.4f}) / {denominator:.4f}\n")
    
    print("--- Final Answer ---")
    print(f"The magnitude of the longitudinal component of the magnetic field Bz is: {abs(Bz):.4e} T")
    
    # Return the final value for the answer block
    return abs(Bz)

# Execute the calculation
final_answer = calculate_magnetic_field_component()