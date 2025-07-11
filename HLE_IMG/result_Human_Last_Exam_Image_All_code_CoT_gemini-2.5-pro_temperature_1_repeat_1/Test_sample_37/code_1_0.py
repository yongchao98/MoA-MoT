import numpy as np

def calculate_bz_component():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element on a cylinder surface, based on the Biot-Savart law.
    """
    # --- Given Parameters ---
    # Permeability of free space (H/m or T*m/A)
    mu0 = 4 * np.pi * 1e-7
    
    # Source point in cylindrical coordinates
    r_s_cyl = 0.34  # m
    theta_s_deg = 60  # degrees
    z_s_cyl = 0.5   # m
    
    # Target field point in Cartesian coordinates
    x_f = 0.1   # m
    y_f = -0.1  # m
    z_f = 0.2   # m
    
    # Azimuthal surface current density magnitude.
    # Assuming the given "J" in A/m^2 is a typo for surface current K in A/m.
    K_mag = 40  # A/m

    # --- Calculations ---
    
    # Step 1: Convert source point to Cartesian coordinates
    theta_s_rad = np.deg2rad(theta_s_deg)
    x_s = r_s_cyl * np.cos(theta_s_rad)
    y_s = r_s_cyl * np.sin(theta_s_rad)
    z_s = z_s_cyl
    r_s_vec = np.array([x_s, y_s, z_s])
    
    # Step 2: Define target field point vector
    r_f_vec = np.array([x_f, y_f, z_f])
    
    # Step 3: Calculate the displacement vector R and its magnitude
    R_vec = r_f_vec - r_s_vec
    R_mag = np.linalg.norm(R_vec)
    
    # Step 4: Define the source current vector K in Cartesian coordinates
    # The current is azimuthal, so its direction is theta_hat.
    # In Cartesian coordinates, theta_hat = (-sin(theta), cos(theta), 0).
    K_x = K_mag * (-np.sin(theta_s_rad))
    K_y = K_mag * (np.cos(theta_s_rad))
    K_z = 0
    K_vec = np.array([K_x, K_y, K_z])
    
    # Step 5: Calculate the cross product K x R
    cross_product_vec = np.cross(K_vec, R_vec)
    cross_product_z = cross_product_vec[2]
    
    # Step 6: Apply the Biot-Savart law to find the Bz value
    # We are calculating the integrand value, as the dS area element is not given.
    # Bz = (mu0 / 4pi) * (K x R)_z / R^3
    Bz = (mu0 / (4 * np.pi)) * cross_product_z / (R_mag**3)

    # --- Print the detailed steps and the final equation ---
    print("This calculation is based on the Biot-Savart law: B_z = (μ₀ / 4π) * ([K × R]_z / R³)")
    print("Assuming the given current density J=40 A/m² is a surface current density K=40 A/m.\n")

    print("Step 1: Source and Field Points in Cartesian Coordinates")
    print(f"Source point r_s = ({r_s_vec[0]:.4f}, {r_s_vec[1]:.4f}, {r_s_vec[2]:.4f}) m")
    print(f"Field point r_f = ({r_f_vec[0]:.4f}, {r_f_vec[1]:.4f}, {r_f_vec[2]:.4f}) m\n")
    
    print("Step 2: Displacement Vector R")
    print(f"R = r_f - r_s = ({R_vec[0]:.4f}, {R_vec[1]:.4f}, {R_vec[2]:.4f}) m")
    print(f"Magnitude |R| = {R_mag:.4f} m\n")

    print("Step 3: Current Density Vector K")
    print(f"K = ({K_vec[0]:.4f}, {K_vec[1]:.4f}, {K_vec[2]:.4f}) A/m\n")

    print("Step 4: z-Component of the Cross Product [K × R]")
    print(f"[K × R]_z = K_x * R_y - K_y * R_x")
    print(f"[K × R]_z = ({K_vec[0]:.4f}) * ({R_vec[1]:.4f}) - ({K_vec[1]:.4f}) * ({R_vec[0]:.4f})")
    print(f"[K × R]_z = {cross_product_z:.4f} A·m\n")

    print("Step 5: Final Calculation for B_z")
    print(f"The term μ₀ / 4π is equal to {1e-7:.1e} T·m/A.")
    print("The final equation for B_z is:")
    print(f"B_z = {1e-7:.1e} * ({cross_product_z:.4f} / {R_mag:.4f}³)")
    print(f"B_z = {1e-7:.1e} * ({cross_product_z:.4f} / {R_mag**3:.4f})")
    print(f"B_z = {Bz:.4e} T")
    
    return Bz

if __name__ == '__main__':
    final_answer = calculate_bz_component()
    # The final answer is returned in scientific notation.
    # For example, 1.2017e-05 is 0.000012017
    # print(f"\nFinal Answer: {final_answer}")
    # <<<1.2017e-05>>> # This format is for the final answer extraction