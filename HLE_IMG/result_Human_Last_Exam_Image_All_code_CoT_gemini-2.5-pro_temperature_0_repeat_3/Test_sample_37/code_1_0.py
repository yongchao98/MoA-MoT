import numpy as np

def calculate_magnetic_field():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current density at a source point on a cylinder surface.
    """
    # --- Given Parameters ---
    mu_0 = 4 * np.pi * 1e-7  # Permeability of free space (H/m)
    J_theta_mag = 40.0       # Magnitude of azimuthal current density (A/m^2)

    # Target field point in Cartesian coordinates (m)
    r_f = np.array([0.1, -0.1, 0.2])

    # Source point in Cylindrical coordinates (m, degrees, m)
    r_s_cyl, theta_s_deg, z_s_cyl = 0.34, 60.0, 0.5

    # --- Step 1: Convert source point to Cartesian coordinates ---
    theta_s_rad = np.radians(theta_s_deg)
    xs = r_s_cyl * np.cos(theta_s_rad)
    ys = r_s_cyl * np.sin(theta_s_rad)
    zs = z_s_cyl
    r_s = np.array([xs, ys, zs])

    # --- Step 2: Calculate the displacement vector R and its magnitude ---
    R = r_f - r_s
    R_mag = np.linalg.norm(R)

    # --- Step 3: Determine the current density vector J in Cartesian coordinates ---
    # The direction is azimuthal (a_theta). In Cartesian coordinates:
    # a_theta = -sin(theta) * i + cos(theta) * j
    a_theta = np.array([-np.sin(theta_s_rad), np.cos(theta_s_rad), 0])
    J = J_theta_mag * a_theta

    # --- Step 4: Calculate the cross product J x R ---
    cross_product = np.cross(J, R)
    cross_product_z = cross_product[2]

    # --- Step 5: Calculate the longitudinal component of the magnetic field (Bz) ---
    # Using the Biot-Savart law for a volume current element.
    # We assume a unit volume element (dV = 1 m^3) to find the field in Tesla.
    # B = (μ₀ / 4π) * (J x R) / |R|^3
    Bz = (mu_0 / (4 * np.pi)) * cross_product_z / (R_mag**3)

    # --- Step 6: Output the final equation with all numbers ---
    print("The final equation for the longitudinal component of the magnetic field (B_z) is:")
    print("B_z = (μ₀ / 4π) * [J x R]_z / |R|³")
    print(f"B_z = (1e-7) * ({cross_product_z:.4f}) / ({R_mag:.4f})³")
    print(f"B_z = (1e-7) * ({cross_product_z:.4f}) / ({R_mag**3:.4f})")
    print(f"B_z = {Bz:.4e} T")

if __name__ == '__main__':
    calculate_magnetic_field()