import math

def calculate_bz():
    """
    Calculates the longitudinal component of the magnetic field (Bz) from a source
    point on a cylinder to a target field point.
    """
    # --- Given Parameters ---
    # Permeability of free space constant
    mu0_over_4pi = 1e-7  # H/m

    # Cylinder parameters (not all are needed for this specific calculation)
    # r_cyl = 0.34  # m
    # L_cyl = 1.4   # m

    # Source point in cylindrical coordinates (r, theta, z)
    r_s_cyl = 0.34       # m
    theta_s_deg = 60.0   # degrees
    z_s = 0.5          # m

    # Target field point in Cartesian coordinates (x, y, z)
    x_f = 0.1          # m
    y_f = -0.1         # m
    z_f = 0.2          # m

    # Azimuthal current density.
    # The unit A/m^2 is inconsistent with surface current density (A/m).
    # Assuming it's the magnitude of surface current density K_theta.
    K_theta = 40.0     # A/m

    # The problem asks for the field B (in Tesla) from a point source.
    # This implies a discrete element calculation. Assuming a source area of 1 m^2.
    delta_S = 1.0      # m^2

    # --- Step 1: Convert source point to Cartesian coordinates ---
    theta_s_rad = math.radians(theta_s_deg)
    x_s = r_s_cyl * math.cos(theta_s_rad)
    y_s = r_s_cyl * math.sin(theta_s_rad)
    # z_s is already given

    # --- Step 2: Calculate the displacement vector R = r_f - r_s ---
    Rx = x_f - x_s
    Ry = y_f - y_s
    Rz = z_f - z_s

    # --- Step 3: Calculate the magnitude of R ---
    R_mag_sq = Rx**2 + Ry**2 + Rz**2
    R_mag = math.sqrt(R_mag_sq)

    # --- Step 4: Convert current density vector K to Cartesian coordinates ---
    # K = K_theta * (-sin(theta) i_hat + cos(theta) j_hat)
    Kx = -K_theta * math.sin(theta_s_rad)
    Ky = K_theta * math.cos(theta_s_rad)
    # Kz = 0

    # --- Step 5: Calculate the z-component of the cross product (K x R) ---
    cross_product_z = Kx * Ry - Ky * Rx

    # --- Step 6: Calculate Bz ---
    Bz = mu0_over_4pi * cross_product_z / (R_mag**3) * delta_S
    
    # --- Step 7: Print the final equation with numerical values ---
    print("The longitudinal component of the magnetic field (Bz) is calculated using the Biot-Savart law for a discrete current element:")
    print("Bz = (μ₀/4π) * (K x R)_z / |R|³ * ΔS")
    print("\nSubstituting the numerical values:")
    print(f"Bz = ({mu0_over_4pi:.1e}) * (({Kx:.4f}) * ({Ry:.4f}) - ({Ky:.4f}) * ({Rx:.4f})) / ({R_mag:.4f})³ * {delta_S}")
    print(f"Bz = ({mu0_over_4pi:.1e}) * ({cross_product_z:.4f}) / ({R_mag**3:.4f}) * {delta_S}")
    print(f"Bz = {Bz:.4e} T")

if __name__ == '__main__':
    calculate_bz()
<<<1.2019e-05>>>