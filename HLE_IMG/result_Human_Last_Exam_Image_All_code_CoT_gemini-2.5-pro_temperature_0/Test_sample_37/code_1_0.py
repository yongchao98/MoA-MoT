import math

def calculate_bz():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current density at a source point on a cylinder.
    """
    # Step 1: Define constants and given variables
    mu_0 = 4 * math.pi * 1e-7  # Permeability of free space (H/m)
    mu_0_over_4pi = 1e-7

    # Target field point in Cartesian coordinates (m)
    x_f, y_f, z_f = 0.1, -0.1, 0.2
    r_f = (x_f, y_f, z_f)

    # Source point in Cylindrical coordinates
    r_s_cyl, theta_s_deg, z_s_cyl = 0.34, 60, 0.5 # (m, degrees, m)

    # Azimuthal volume current density (A/m^2)
    J_theta_mag = 40.0

    print("Given Parameters:")
    print(f"Target field point (x_f, y_f, z_f) = {r_f} m")
    print(f"Source point (r, θ, z) = ({r_s_cyl} m, {theta_s_deg}°, {z_s_cyl} m)")
    print(f"Azimuthal current density magnitude |J_θ| = {J_theta_mag} A/m^2")
    print(f"Permeability of free space / 4π (μ₀/4π) = {mu_0_over_4pi} H/m\n")

    # Step 2: Convert source point to Cartesian coordinates
    theta_s_rad = math.radians(theta_s_deg)
    x_s = r_s_cyl * math.cos(theta_s_rad)
    y_s = r_s_cyl * math.sin(theta_s_rad)
    z_s = z_s_cyl
    r_s_cart = (x_s, y_s, z_s)
    
    # Step 3: Calculate the displacement vector R
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    R_vec = (R_x, R_y, R_z)
    R_mag = math.sqrt(R_x**2 + R_y**2 + R_z**2)

    # Step 4: Convert the current density vector J to Cartesian coordinates
    # The azimuthal unit vector θ_hat = -sin(θ)i + cos(θ)j
    J_x = -J_theta_mag * math.sin(theta_s_rad)
    J_y = J_theta_mag * math.cos(theta_s_rad)
    J_z = 0
    J_vec = (J_x, J_y, J_z)

    # Step 5 & 6: Calculate the z-component of the cross product (J x R)
    cross_product_z = J_x * R_y - J_y * R_x

    # Step 7: Calculate Bz using the Biot-Savart Law
    # B_z = (μ₀/4π) * (J × R)_z / |R|³
    R_mag_cubed = R_mag**3
    B_z = mu_0_over_4pi * cross_product_z / R_mag_cubed

    print("Calculation Steps:")
    print(f"1. Source point in Cartesian (x_s, y_s, z_s) = ({x_s:.4f}, {y_s:.4f}, {z_s:.4f}) m")
    print(f"2. Displacement vector R = r_f - r_s = ({R_x:.4f}, {R_y:.4f}, {R_z:.4f}) m")
    print(f"3. Magnitude of displacement vector |R| = {R_mag:.4f} m")
    print(f"4. Current density vector J = ({J_x:.4f}, {J_y:.4f}, {J_z:.4f}) A/m^2")
    print(f"5. z-component of cross product (J × R)_z = {cross_product_z:.4f} A/m")
    
    print("\nFinal Equation for B_z:")
    print(f"B_z = (μ₀/4π) * (J × R)_z / |R|³")
    print(f"B_z = ({mu_0_over_4pi}) * ({cross_product_z:.4f}) / ({R_mag:.4f})³")
    print(f"B_z = ({mu_0_over_4pi}) * ({cross_product_z:.4f}) / ({R_mag_cubed:.4f})")
    print(f"B_z = {B_z:.4e} T")
    
    return B_z

# Run the calculation and print the final answer
final_B_z = calculate_bz()
# The final answer is requested in a specific format.
# The magnitude is the absolute value of the calculated field.
final_answer = abs(final_B_z)
print(f"\nThe magnitude of the longitudinal component of the magnetic field B_z is {final_answer:.4e} T.")
print(f"<<<{final_answer:.4e}>>>")
