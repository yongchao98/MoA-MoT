import numpy as np

def calculate_bz_component():
    """
    Calculates the longitudinal component of the magnetic field (B_z) at a target point
    due to a current density at a source point on a cylinder.
    """
    # --- Given Constants and Values ---
    # Permeability of free space / 4*pi (T*m/A)
    mu0_over_4pi = 1e-7
    
    # Source point in cylindrical coordinates (r, theta_deg, z)
    r_s = 0.34  # meters
    theta_deg_s = 60  # degrees
    z_s = 0.5  # meters
    
    # Target field point in Cartesian coordinates (x_f, y_f, z_f)
    xf, yf, zf = 0.1, -0.1, 0.2
    rf_vec = np.array([xf, yf, zf])
    
    # Azimuthal surface current density K_theta (A/m)
    # Assuming the given "40 A/m^2" is a typo for surface current density K=40 A/m.
    K_theta = 40.0
    
    # Assumed differential surface area (m^2) to get a finite answer in Tesla.
    dS = 1.0
    
    # --- Step 1: Convert source point to Cartesian coordinates ---
    theta_rad_s = np.deg2rad(theta_deg_s)
    xs = r_s * np.cos(theta_rad_s)
    ys = r_s * np.sin(theta_rad_s)
    rs_vec = np.array([xs, ys, z_s])
    
    # --- Step 2: Calculate the displacement vector R ---
    R_vec = rf_vec - rs_vec
    R_mag = np.linalg.norm(R_vec)
    
    # --- Step 3: Define the current density vector K in Cartesian coordinates ---
    # The azimuthal unit vector theta_hat = [-sin(theta), cos(theta), 0]
    theta_hat = np.array([-np.sin(theta_rad_s), np.cos(theta_rad_s), 0])
    K_vec = K_theta * theta_hat
    
    # --- Step 4: Calculate the cross product (K x R) ---
    K_cross_R = np.cross(K_vec, R_vec)
    K_cross_R_z = K_cross_R[2]
    
    # --- Step 5: Calculate the final B_z value ---
    # B_z = (mu0 / 4pi) * (K_cross_R_z / R^3) * dS
    R_mag_cubed = R_mag**3
    Bz = mu0_over_4pi * K_cross_R_z / R_mag_cubed * dS

    # --- Print the results ---
    print("This calculation uses the Biot-Savart Law to find the magnetic field from a current element.")
    print("It assumes the provided value is a surface current density K = 40 A/m, and calculates the field from a unit surface element dS = 1 m².")
    print("\nThe equation for the longitudinal component of the magnetic field is:")
    print("B_z = (μ₀/4π) * (K × R)z / |R|³ * dS")
    print("\nPlugging in the calculated values:")
    print(f"B_z = ({mu0_over_4pi:.1e}) * ({K_cross_R_z:.4f}) / ({R_mag_cubed:.4f}) * {dS}")
    print(f"\nThe calculated magnitude of the longitudinal component of the magnetic field is:")
    print(f"B_z = {Bz:.4e} T")

# Run the calculation
calculate_bz_component()