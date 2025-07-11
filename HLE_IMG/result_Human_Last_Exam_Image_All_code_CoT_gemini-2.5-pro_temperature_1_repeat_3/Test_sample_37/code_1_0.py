import math

def calculate_magnetic_field_component():
    """
    Calculates the longitudinal component of the magnetic field (B_z)
    at a target point due to a current density on a cylinder surface.
    """
    # 1. Define Constants and Given Values
    mu_0 = 4 * math.pi * 1e-7  # Permeability of free space in H/m
    mu_0_over_4pi = 1e-7

    # Target field point in Cartesian coordinates (m)
    x_f, y_f, z_f = 0.1, -0.1, 0.2

    # Source point in Cylindrical coordinates (m, degrees)
    r_s_cyl, theta_s_deg, z_s = 0.34, 60.0, 0.5

    # Azimuthal current density (A/m^2)
    J_theta = 40.0

    # 2. Convert Source Point and Current Density to Cartesian Coordinates
    theta_s_rad = math.radians(theta_s_deg)

    # Source point in Cartesian (x_s, y_s, z_s)
    x_s = r_s_cyl * math.cos(theta_s_rad)
    y_s = r_s_cyl * math.sin(theta_s_rad)

    # Current density vector in Cartesian (J_x, J_y, J_z)
    # The azimuthal unit vector theta_hat is (-sin(theta), cos(theta), 0)
    J_x = -J_theta * math.sin(theta_s_rad)
    J_y = J_theta * math.cos(theta_s_rad)
    J_z = 0.0

    # 3. Calculate the Displacement Vector R and its Magnitude
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    R_mag = math.sqrt(R_x**2 + R_y**2 + R_z**2)

    # 4. Calculate the z-component of the cross product (J x R)
    J_cross_R_z = J_x * R_y - J_y * R_x

    # 5. Calculate the final B_z value
    # The formula is dB_z = (mu_0 / 4pi) * (J dS x R)_z / |R|^3.
    # We assume dS = 1 m^2 to resolve the units.
    B_z = mu_0_over_4pi * J_cross_R_z / (R_mag**3)

    # --- Print the detailed calculation ---
    print("--- Calculation of the Magnetic Field Component B_z ---")
    print("\nThe Biot-Savart Law for a surface current element is B = (mu_0 / 4pi) * (J dS x R) / |R|^3.")
    print("Since J is given in A/m^2, we assume a unit area dS = 1 m^2 to find the field contribution in Tesla.")

    print("\nStep 1: Determine source and displacement vectors in Cartesian coordinates.")
    print(f"Source point (r={r_s_cyl} m, theta={theta_s_deg} deg, z={z_s} m) is ({x_s:.5f}, {y_s:.5f}, {z_s:.1f}) m in Cartesian.")
    print(f"Target point is ({x_f}, {y_f}, {z_f}) m.")
    print(f"Displacement vector R = Target - Source = ({R_x:.5f}, {R_y:.5f}, {R_z:.5f}) m.")
    print(f"Magnitude |R| = {R_mag:.5f} m.")

    print("\nStep 2: Determine the current density vector J in Cartesian coordinates.")
    print(f"Current density J (azimuthal) = ({J_x:.5f}, {J_y:.5f}, {J_z:.1f}) A/m^2.")

    print("\nStep 3: Calculate the z-component of the cross product (J x R).")
    print(f"(J x R)_z = J_x * R_y - J_y * R_x")
    print(f"(J x R)_z = ({J_x:.5f}) * ({R_y:.5f}) - ({J_y:.5f}) * ({R_x:.5f})")
    print(f"(J x R)_z = {J_cross_R_z:.5f} A/m.")

    print("\nStep 4: Calculate B_z using the Biot-Savart Law.")
    print(f"B_z = (mu_0 / 4pi) * (J x R)_z / |R|^3")
    print(f"B_z = ({mu_0_over_4pi:.2e}) * ({J_cross_R_z:.5f}) / ({R_mag:.5f})^3")
    print(f"B_z = ({mu_0_over_4pi:.2e}) * ({J_cross_R_z:.5f}) / ({R_mag**3:.5f})")
    print(f"B_z = {B_z:.4e} T")

    print(f"\nThe magnitude of the longitudinal component of the magnetic field (B_z) is {B_z:.4e} T.")
    print(f"<<<{B_z:.4e}>>>")

calculate_magnetic_field_component()