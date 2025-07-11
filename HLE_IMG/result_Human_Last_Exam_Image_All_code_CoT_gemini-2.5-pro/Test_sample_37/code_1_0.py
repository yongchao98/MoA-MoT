import numpy as np

def calculate_magnetic_field_component():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element on the surface of a cylinder.
    """
    # Step 1: Define parameters and constants in SI units
    # Permeability of free space
    mu_0 = 4 * np.pi * 1e-7  # H/m
    # The term mu_0 / (4 * pi) simplifies to 1e-7
    mu_0_over_4pi = 1e-7

    # Current Density: Assuming the provided "J" is a surface current density "K".
    # The unit A/m^2 is assumed to be a typo for A/m.
    K_theta_mag = 40.0  # A/m

    # To calculate a finite B-field value from a point source contribution,
    # we assume a unit surface area element.
    dS = 1.0  # m^2

    # Target field point in Cartesian coordinates
    xf, yf, zf = 0.1, -0.1, 0.2
    r_f = np.array([xf, yf, zf])

    # Source point in Cylindrical coordinates
    r_s_cyl = 0.34
    theta_s_deg = 60.0
    z_s_cyl = 0.5
    # Convert angle to radians for trigonometric functions
    theta_s_rad = np.deg2rad(theta_s_deg)

    # Step 2: Convert source point to Cartesian coordinates
    xs = r_s_cyl * np.cos(theta_s_rad)
    ys = r_s_cyl * np.sin(theta_s_rad)
    zs = z_s_cyl
    r_s = np.array([xs, ys, zs])

    # Step 3: Calculate the displacement vector R and its magnitude
    R_vec = r_f - r_s
    R_mag = np.linalg.norm(R_vec)

    # Step 4: Express the azimuthal current density vector K in Cartesian coordinates
    # The azimuthal unit vector (theta_hat) in Cartesian coordinates is (-sin(theta), cos(theta), 0)
    K_vec_cartesian = K_theta_mag * np.array([-np.sin(theta_s_rad), np.cos(theta_s_rad), 0])

    # Step 5: Calculate the cross product K x R
    K_cross_R = np.cross(K_vec_cartesian, R_vec)
    K_cross_R_z = K_cross_R[2]

    # Step 6: Apply the Biot-Savart law to find the z-component of the magnetic field
    # B_z = (mu_0 / 4π) * (K x R)_z / |R|^3 * dS
    B_z = mu_0_over_4pi * K_cross_R_z / (R_mag**3) * dS

    # --- Print the detailed calculation ---
    print("--- Calculation of the Longitudinal Magnetic Field (Bz) ---")
    print("\nStep 1: Define problem parameters")
    print(f"Target point r_f = <{xf}, {yf}, {zf}> m")
    print(f"Source point (cylindrical) = ({r_s_cyl} m, {theta_s_deg}°, {z_s_cyl} m)")
    print(f"Surface current density K_θ = {K_theta_mag} A/m")

    print("\nStep 2: Convert source point to Cartesian coordinates r_s")
    print(f"r_s = <{xs:.4f}, {ys:.4f}, {zs:.4f}> m")

    print("\nStep 3: Calculate displacement vector R = r_f - r_s")
    print(f"R = <{R_vec[0]:.4f}, {R_vec[1]:.4f}, {R_vec[2]:.4f}> m")
    print(f"|R| = {R_mag:.4f} m")

    print("\nStep 4: Express current density vector K in Cartesian coordinates")
    print(f"K = <{K_vec_cartesian[0]:.4f}, {K_vec_cartesian[1]:.4f}, {K_vec_cartesian[2]:.4f}> A/m")

    print("\nStep 5: Calculate the final equation for Bz")
    print("B_z = (μ₀ / 4π) * (K x R)_z / |R|³")
    print(f"B_z = ({mu_0_over_4pi:.1e}) * ({K_cross_R_z:.4f}) / ({R_mag:.4f})³")
    print(f"B_z = ({mu_0_over_4pi:.1e}) * ({K_cross_R_z:.4f}) / ({R_mag**3:.4f})")
    print(f"B_z = {B_z:.4e} T")
    print("\nNote: The result is for a representative surface element dS = 1 m².")

calculate_magnetic_field_component()
<<<1.2017e-05>>>