import numpy as np

def calculate_magnetic_field():
    """
    Calculates the longitudinal component of the magnetic field (B_z) at a target point
    due to a current element at a source point on a cylinder.
    """
    # --- Given Parameters ---
    # Permeability of free space
    mu_0 = 4 * np.pi * 1e-7  # H/m

    # Source point in cylindrical coordinates
    r_s_cyl = 0.34  # m
    theta_s_deg = 60.0  # degrees
    z_s = 0.5  # m

    # Target field point in Cartesian coordinates
    x_f = 0.1  # m
    y_f = -0.1  # m
    z_f = 0.2  # m

    # Azimuthal current density (interpreted as volume density J due to A/m^2 units)
    J_theta = 40.0  # A/m^2

    # --- Calculations ---

    # 1. Convert source point to Cartesian coordinates
    theta_s_rad = np.deg2rad(theta_s_deg)
    x_s = r_s_cyl * np.cos(theta_s_rad)
    y_s = r_s_cyl * np.sin(theta_s_rad)

    # 2. Represent the current density vector J in Cartesian coordinates
    # The azimuthal unit vector theta_hat in Cartesian is (-sin(theta), cos(theta), 0)
    J_x = -J_theta * np.sin(theta_s_rad)
    J_y = J_theta * np.cos(theta_s_rad)
    J_z = 0.0

    # 3. Calculate the displacement vector R from source to target
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s

    # 4. Calculate the magnitude of R
    R_mag_sq = R_x**2 + R_y**2 + R_z**2
    R_mag = np.sqrt(R_mag_sq)

    # 5. Calculate the z-component of the cross product (J x R)
    cross_product_z = J_x * R_y - J_y * R_x

    # 6. Calculate B_z using the Biot-Savart Law for a unit volume element (dV=1)
    # The factor (mu_0 / 4pi) simplifies to 1e-7
    B_z = (mu_0 / (4 * np.pi)) * (cross_product_z / (R_mag**3))

    # --- Output Results ---
    print("--- Input Parameters ---")
    print(f"Source Point (cylindrical): r={r_s_cyl} m, theta={theta_s_deg} deg, z={z_s} m")
    print(f"Target Point (Cartesian): x={x_f} m, y={y_f} m, z={z_f} m")
    print(f"Current Density (azimuthal): J_theta = {J_theta} A/m^2")
    print(f"Permeability of Free Space (mu_0): {mu_0:.2e} H/m\n")

    print("--- Calculation Steps ---")
    print("1. Source Point in Cartesian Coordinates:")
    print(f"   x_s = {r_s_cyl} * cos({theta_s_deg}°) = {x_s:.4f} m")
    print(f"   y_s = {r_s_cyl} * sin({theta_s_deg}°) = {y_s:.4f} m\n")

    print("2. Current Density in Cartesian Coordinates:")
    print(f"   J_x = -{J_theta} * sin({theta_s_deg}°) = {J_x:.4f} A/m^2")
    print(f"   J_y = {J_theta} * cos({theta_s_deg}°) = {J_y:.4f} A/m^2\n")

    print("3. Displacement Vector (R = r_f - r_s):")
    print(f"   R_x = {x_f} - {x_s:.4f} = {R_x:.4f} m")
    print(f"   R_y = {y_f} - {y_s:.4f} = {R_y:.4f} m")
    print(f"   R_z = {z_f} - {z_s:.4f} = {R_z:.4f} m\n")

    print("4. Magnitude of Displacement Vector (R):")
    print(f"   R = sqrt(({R_x:.4f})^2 + ({R_y:.4f})^2 + ({R_z:.4f})^2) = {R_mag:.4f} m\n")

    print("5. Z-Component of Cross Product (J x R)_z:")
    print(f"   (J x R)_z = J_x*R_y - J_y*R_x = ({J_x:.4f})*({R_y:.4f}) - ({J_y:.4f})*({R_x:.4f}) = {cross_product_z:.4f} A/m\n")

    print("--- Final Equation for B_z ---")
    print("B_z = (mu_0 / 4pi) * ( (J_x * R_y - J_y * R_x) / (R_x^2 + R_y^2 + R_z^2)^(3/2) )")
    print(f"B_z = (1e-7) * ( ({J_x:.4f} * {R_y:.4f} - {J_y:.4f} * {R_x:.4f}) / (({R_x:.4f})^2 + ({R_y:.4f})^2 + ({R_z:.4f})^2)^(3/2) )")
    print(f"B_z = (1e-7) * ( {cross_product_z:.4f} / ({R_mag:.4f})^3 )")
    print(f"B_z = (1e-7) * ( {cross_product_z:.4f} / {R_mag**3:.4f} )")
    print(f"B_z = {B_z:.4e} T\n")

    print(f"The magnitude of the longitudinal component of the magnetic field B_z is {B_z:.4e} Tesla.")


if __name__ == '__main__':
    calculate_magnetic_field()
    # To get the final value for the submission format
    mu_0 = 4 * np.pi * 1e-7
    r_s_cyl = 0.34
    theta_s_deg = 60.0
    z_s = 0.5
    x_f = 0.1
    y_f = -0.1
    z_f = 0.2
    J_theta = 40.0
    theta_s_rad = np.deg2rad(theta_s_deg)
    x_s = r_s_cyl * np.cos(theta_s_rad)
    y_s = r_s_cyl * np.sin(theta_s_rad)
    J_x = -J_theta * np.sin(theta_s_rad)
    J_y = J_theta * np.cos(theta_s_rad)
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    R_mag = np.sqrt(R_x**2 + R_y**2 + R_z**2)
    cross_product_z = J_x * R_y - J_y * R_x
    B_z = (mu_0 / (4 * np.pi)) * (cross_product_z / (R_mag**3))
    print(f'<<<{B_z:.3e}>>>')