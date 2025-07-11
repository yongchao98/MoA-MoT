import numpy as np

def calculate_bz():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element on the surface of a cylinder.
    """
    # --- 1. Define constants and given values ---
    mu_0_over_4pi = 1e-7  # Permeability of free space constant (T*m/A)
    
    # Cylinder dimensions (not directly needed for this point-source calculation but given for context)
    r_cylinder = 0.34  # meters
    L_cylinder = 1.4   # meters
    
    # Target field point (Cartesian coordinates)
    xf, yf, zf = 0.1, -0.1, 0.2
    r_f = np.array([xf, yf, zf])
    
    # Source point (Cylindrical coordinates)
    r_s_cyl, theta_s_deg, z_s = 0.34, 60.0, 0.5
    
    # Azimuthal surface current density (assuming A/m, typo in problem's A/m^2)
    K_theta = 40.0  # A/m
    
    print("Step 1: Convert source point from cylindrical to Cartesian coordinates.")
    # Convert theta from degrees to radians
    theta_s_rad = np.deg2rad(theta_s_deg)
    
    # Calculate Cartesian coordinates of the source point
    x_s = r_s_cyl * np.cos(theta_s_rad)
    y_s = r_s_cyl * np.sin(theta_s_rad)
    # z_s is already given
    r_s = np.array([x_s, y_s, z_s])
    print(f"Source point (r, θ, z) = ({r_s_cyl} m, {theta_s_deg}°, {z_s} m)")
    print(f"Source point (x, y, z) = ({x_s:.4f}, {y_s:.4f}, {z_s:.4f}) m\n")

    print("Step 2: Calculate the displacement vector R from source to target.")
    # Calculate displacement vector R = r_f - r_s
    R_vec = r_f - r_s
    print(f"Target vector r_f = {r_f} m")
    print(f"Source vector r_s = [{r_s[0]:.4f}, {r_s[1]:.4f}, {r_s[2]:.4f}] m")
    print(f"Displacement vector R = r_f - r_s = [{R_vec[0]:.4f}, {R_vec[1]:.4f}, {R_vec[2]:.4f}] m\n")

    print("Step 3: Convert the current density vector K from cylindrical to Cartesian coordinates.")
    # K is purely azimuthal (in the theta-hat direction)
    # The Cartesian representation of theta-hat is (-sin(theta), cos(theta), 0)
    Kx = -K_theta * np.sin(theta_s_rad)
    Ky = K_theta * np.cos(theta_s_rad)
    Kz = 0.0
    K_vec = np.array([Kx, Ky, Kz])
    print(f"Surface current density K_theta = {K_theta} A/m")
    print(f"Current vector K = [{K_vec[0]:.4f}, {K_vec[1]:.4f}, {K_vec[2]:.4f}] A/m\n")

    print("Step 4: Calculate the magnetic field component Bz using the Biot-Savart Law.")
    # Biot-Savart Law for surface current: B = (μ₀/4π) * integral( (K x R) / |R|³ dS )
    # We calculate the value of the integrand for the given points.
    
    # Calculate the z-component of the cross product (K x R)z = Kx*Ry - Ky*Rx
    cross_product_z = K_vec[0] * R_vec[1] - K_vec[1] * R_vec[0]
    
    # Calculate the magnitude of R, then cube it
    R_mag = np.linalg.norm(R_vec)
    R_mag_cubed = R_mag**3
    
    # Calculate Bz
    Bz = mu_0_over_4pi * cross_product_z / R_mag_cubed

    print("The formula for Bz is: Bz = (μ₀/4π) * (Kx * Ry - Ky * Rx) / |R|³")
    print(f"Plugging in the numbers:")
    print(f"Bz = ({mu_0_over_4pi}) * (({K_vec[0]:.4f}) * ({R_vec[1]:.4f}) - ({K_vec[1]:.4f}) * ({R_vec[0]:.4f})) / ({R_mag:.4f})³")
    print(f"Bz = ({mu_0_over_4pi}) * ({cross_product_z:.4f}) / ({R_mag_cubed:.4f})")
    print(f"Bz = {mu_0_over_4pi * cross_product_z / R_mag_cubed:.4e} T\n")

    print(f"The final calculated magnitude of the longitudinal component of the magnetic field (Bz) is {Bz:.4e} T.")

calculate_bz()