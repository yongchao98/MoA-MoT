import numpy as np

def solve_magnetic_field():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element on a cylinder's surface using the Biot-Savart Law.
    """
    # --- Step 1: Define constants and given parameters ---
    mu_0 = 4 * np.pi * 1e-7  # Permeability of free space in H/m
    mu_0_over_4pi = 1e-7

    # Source point in cylindrical coordinates (r, theta, z)
    r_s = 0.34  # meters
    theta_s_deg = 60.0  # degrees
    z_s = 0.5  # meters

    # Target field point in Cartesian coordinates (x_f, y_f, z_f)
    xf, yf, zf = 0.1, -0.1, 0.2
    r_f_cart = np.array([xf, yf, zf])

    # Interpreted magnitude of the current element |I*dl| from the given J_theta
    I_dl_magnitude = 40.0  # A.m

    print("--- Input Parameters ---")
    print(f"Permeability constant (mu_0 / 4*pi): {mu_0_over_4pi:.1e} H/m")
    print(f"Source point (cylindrical): (r={r_s} m, theta={theta_s_deg} deg, z={z_s} m)")
    print(f"Target field point (Cartesian): (x={xf} m, y={yf} m, z={zf} m)")
    print(f"Magnitude of the current element |I*dl|: {I_dl_magnitude} A.m\n")

    # --- Step 2: Convert source point to Cartesian coordinates ---
    theta_s_rad = np.deg2rad(theta_s_deg)
    xs = r_s * np.cos(theta_s_rad)
    ys = r_s * np.sin(theta_s_rad)
    zs = z_s
    r_s_cart = np.array([xs, ys, zs])

    print("--- Calculation Steps ---")
    print(f"1. Source point in Cartesian coordinates r_s = (xs, ys, zs):")
    print(f"   r_s = ({xs:.5f}, {ys:.5f}, {zs:.5f}) m\n")

    # --- Step 3: Determine the current element vector I*dl ---
    # The direction is azimuthal (theta_hat). In Cartesian coordinates:
    # theta_hat = [-sin(theta), cos(theta), 0]
    theta_hat_cart = np.array([-np.sin(theta_s_rad), np.cos(theta_s_rad), 0])
    I_dl_vector = I_dl_magnitude * theta_hat_cart
    
    print(f"2. Current element vector I*dl:")
    print(f"   I*dl = |I*dl| * [-sin(theta), cos(theta), 0]")
    print(f"   I*dl = ({I_dl_vector[0]:.5f}, {I_dl_vector[1]:.5f}, {I_dl_vector[2]:.5f}) A.m\n")

    # --- Step 4: Calculate the displacement vector R and its magnitude ---
    R_vector = r_f_cart - r_s_cart
    R_magnitude = np.linalg.norm(R_vector)
    
    print(f"3. Displacement vector R = r_f - r_s:")
    print(f"   R = ({R_vector[0]:.5f}, {R_vector[1]:.5f}, {R_vector[2]:.5f}) m")
    print(f"   Magnitude |R| = {R_magnitude:.5f} m\n")

    # --- Step 5: Compute the cross product (I*dl x R) ---
    cross_product = np.cross(I_dl_vector, R_vector)
    
    print(f"4. Cross product C = (I*dl x R):")
    print(f"   C = ({cross_product[0]:.5f}, {cross_product[1]:.5f}, {cross_product[2]:.5f}) A.m^2\n")

    # --- Step 6: Calculate Bz using the Biot-Savart Law ---
    # Bz = (mu_0 / 4*pi) * C_z / |R|^3
    Bz = mu_0_over_4pi * cross_product[2] / (R_magnitude**3)
    
    print("--- Final Calculation for Bz ---")
    print("The longitudinal component of the magnetic field (Bz) is given by:")
    print("Bz = (mu_0 / 4*pi) * (I*dl x R)_z / |R|^3\n")

    print("Substituting the calculated values into the equation:")
    print(f"Bz = ({mu_0_over_4pi:.1e}) * ({cross_product[2]:.5f}) / ({R_magnitude:.5f})^3")
    print(f"Bz = ({mu_0_over_4pi:.1e}) * ({cross_product[2]:.5f}) / ({R_magnitude**3:.5f})")
    print(f"Bz = {Bz:.5e} Tesla\n")

    print(f"The magnitude of the longitudinal component of the magnetic field is {abs(Bz):.4e} T.")
    
    return abs(Bz)

if __name__ == '__main__':
    final_answer = solve_magnetic_field()
    # The final answer is wrapped for parsing.
    # print(f'<<<{final_answer:.4e}>>>')
    
solve_magnetic_field()