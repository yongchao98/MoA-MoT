import numpy as np

def calculate_bz_component():
    """
    Calculates the longitudinal component of the magnetic field (Bz)
    at a target point due to a surface current element on a cylinder.
    """
    # --- Step 1: Define Constants and Coordinates ---
    # Permeability of free space constant
    mu_0_over_4pi = 1e-7  # T*m/A

    # Target field point (Cartesian)
    xf, yf, zf = 0.1, -0.1, 0.2
    r_f = np.array([xf, yf, zf])

    # Source point (Cylindrical)
    r_s_cyl, theta_s_deg, z_s_cyl = 0.34, 60.0, 0.5

    # Azimuthal surface current density magnitude
    J_theta = 40.0  # A/m^2

    # --- Step 2: Coordinate Transformation for the Source Point ---
    theta_s_rad = np.deg2rad(theta_s_deg)
    xs = r_s_cyl * np.cos(theta_s_rad)
    ys = r_s_cyl * np.sin(theta_s_rad)
    zs = z_s_cyl
    r_s = np.array([xs, ys, zs])

    # --- Step 3: Calculate Displacement Vector R ---
    R_vec = r_f - r_s
    Rx, Ry, Rz = R_vec[0], R_vec[1], R_vec[2]
    R_mag = np.linalg.norm(R_vec)
    R_mag_cubed = R_mag**3

    # --- Step 4: Express Current Density in Cartesian Coordinates ---
    # J = J_theta * hat_theta = J_theta * (-sin(theta) i + cos(theta) j)
    Jx = -J_theta * np.sin(theta_s_rad)
    Jy = J_theta * np.cos(theta_s_rad)

    # --- Step 5 & 6: Apply Biot-Savart Law and Calculate Bz ---
    # The z-component of the cross product (J x R) is (Jx * Ry - Jy * Rx)
    cross_product_z = Jx * Ry - Jy * Rx

    # Calculate Bz using the Biot-Savart law integrand.
    # The problem asks for the field B_z from a source J at a point,
    # which implies we calculate the field contribution per unit source area (dS=1).
    Bz = mu_0_over_4pi * cross_product_z / R_mag_cubed

    # --- Step 7: Final Output ---
    print("The formula for the longitudinal component of the magnetic field B_z is:")
    print(f"B_z = (μ₀/4π) * (Jₓ * Rᵧ - Jᵧ * Rₓ) / (Rₓ² + Rᵧ² + R₂²)^(3/2)\n")

    print("--- Breakdown of the values in the equation ---")
    print(f"Constant term μ₀/4π = {mu_0_over_4pi:.1e} T·m/A\n")

    print(f"Target Point (r_f): ({xf}, {yf}, {zf}) m")
    print(f"Source Point (r_s): (r={r_s_cyl}, θ={theta_s_deg}°, z={z_s_cyl}) -> (x={xs:.4f}, y={ys:.4f}, z={zs:.4f}) m\n")

    print("Displacement Vector (R = r_f - r_s):")
    print(f"Rₓ = {xf:.4f} - {xs:.4f} = {Rx:.4f} m")
    print(f"Rᵧ = {yf:.4f} - {ys:.4f} = {Ry:.4f} m")
    print(f"R₂ = {zf:.4f} - {zs:.4f} = {Rz:.4f} m\n")

    print("Current Density Vector in Cartesian (J):")
    print(f"Jₓ = -{J_theta}*sin({theta_s_deg}°) = {Jx:.4f} A/m²")
    print(f"Jᵧ =  {J_theta}*cos({theta_s_deg}°) = {Jy:.4f} A/m²\n")

    print("Numerator: z-component of (J x R)")
    print(f"(Jₓ * Rᵧ - Jᵧ * Rₓ) = ({Jx:.4f}) * ({Ry:.4f}) - ({Jy:.4f}) * ({Rx:.4f})")
    print(f"                  = {(Jx * Ry):.4f} - ({(Jy * Rx):.4f})")
    print(f"                  = {cross_product_z:.4f} A/m\n")

    print("Denominator: Magnitude of R cubed")
    print(f"R³ = (Rₓ² + Rᵧ² + R₂²)^(3/2) = (({Rx:.4f})² + ({Ry:.4f})² + ({Rz:.4f})²)^(3/2)")
    print(f"   = {R_mag_cubed:.4f} m³\n")

    print("--- Final Result ---")
    print(f"B_z = ({mu_0_over_4pi:.1e}) * ({cross_product_z:.4f}) / ({R_mag_cubed:.4f})")
    print(f"The magnitude of the longitudinal component of the magnetic field is: {np.abs(Bz):.4e} T")

if __name__ == "__main__":
    calculate_bz_component()
