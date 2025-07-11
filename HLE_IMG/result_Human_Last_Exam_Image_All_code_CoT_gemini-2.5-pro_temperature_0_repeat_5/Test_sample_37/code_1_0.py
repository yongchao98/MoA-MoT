import math

def calculate_magnetic_field_component():
    """
    Calculates the longitudinal component of the magnetic field (B_z) at a target point
    due to a current element on a cylinder surface.
    """
    # Step 1: Define constants and given values
    # The term (mu_0 / 4*pi) simplifies to 1e-7
    mu_0_div_4pi = 1e-7  # Units: H/m or T*m/A

    # Target field point in Cartesian coordinates (m)
    x_f, y_f, z_f = 0.1, -0.1, 0.2

    # Source point in Cylindrical coordinates
    r_s_cyl, theta_s_deg, z_s_cyl = 0.34, 60.0, 0.5  # (m, degrees, m)

    # The problem states "azimuthal current density (J_theta) ... is 40 A/m^2".
    # This is dimensionally inconsistent for a surface current. To obtain a magnetic field
    # in Tesla, we interpret this as the magnitude of the current element, |I*dl| = 40 A*m.
    I_dl_magnitude = 40.0  # Units: A*m

    # Step 2: Convert source point to Cartesian coordinates
    # Convert angle from degrees to radians
    theta_s_rad = math.radians(theta_s_deg)

    # Calculate Cartesian coordinates of the source point
    x_s = r_s_cyl * math.cos(theta_s_rad)
    y_s = r_s_cyl * math.sin(theta_s_rad)
    z_s = z_s_cyl

    # Step 3: Calculate the displacement vector R = r_f - r_s
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s

    # Step 4: Calculate the magnitude of R
    R_magnitude = math.sqrt(R_x**2 + R_y**2 + R_z**2)

    # Step 5: Determine the current element vector I_dl in Cartesian coordinates
    # The azimuthal unit vector (theta-hat) in Cartesian coordinates is (-sin(theta), cos(theta), 0)
    I_dl_x = I_dl_magnitude * (-math.sin(theta_s_rad))
    I_dl_y = I_dl_magnitude * (math.cos(theta_s_rad))

    # Step 6: Calculate the z-component of the cross product (I_dl x R)
    cross_product_z = I_dl_x * R_y - I_dl_y * R_x

    # Step 7: Apply the Biot-Savart law for B_z
    # B_z = (mu_0 / 4*pi) * (I_dl x R)_z / |R|^3
    B_z = mu_0_div_4pi * cross_product_z / (R_magnitude**3)

    # Step 8: Print the final equation and result
    print("Calculation of the longitudinal magnetic field component (B_z)")
    print("------------------------------------------------------------")
    print("The Biot-Savart Law for the z-component is: B_z = (μ₀ / 4π) * (I·dl × R)_z / |R|³")
    
    print("\nGiven values:")
    print(f"  Target point (x_f, y_f, z_f) = ({x_f}, {y_f}, {z_f}) m")
    print(f"  Source point (r, θ, z) = ({r_s_cyl} m, {theta_s_deg}°, {z_s_cyl} m)")
    print(f"  Current element magnitude |I·dl| = {I_dl_magnitude} A·m (interpreted from the problem statement)")
    print(f"  μ₀ / 4π = {mu_0_div_4pi:.1e} T·m/A")

    print("\nStep 1: Source point in Cartesian coordinates (x_s, y_s, z_s)")
    print(f"  x_s = r·cos(θ) = {r_s_cyl}·cos({theta_s_deg}°) = {x_s:.4f} m")
    print(f"  y_s = r·sin(θ) = {r_s_cyl}·sin({theta_s_deg}°) = {y_s:.4f} m")
    print(f"  z_s = {z_s:.4f} m")

    print("\nStep 2: Displacement vector R = r_f - r_s")
    print(f"  R_x = x_f - x_s = {x_f} - {x_s:.4f} = {R_x:.4f} m")
    print(f"  R_y = y_f - y_s = {y_f} - {y_s:.4f} = {R_y:.4f} m")
    print(f"  R_z = z_f - z_s = {z_f} - {z_s:.4f} = {R_z:.4f} m")
    print(f"  |R| = sqrt(R_x² + R_y² + R_z²) = sqrt({R_x:.4f}² + {R_y:.4f}² + {R_z:.4f}²) = {R_magnitude:.4f} m")

    print("\nStep 3: Current element vector I·dl in Cartesian coordinates")
    print(f"  I·dl_x = |I·dl|·(-sin(θ)) = {I_dl_magnitude}·(-sin({theta_s_deg}°)) = {I_dl_x:.4f} A·m")
    print(f"  I·dl_y = |I·dl|·(cos(θ)) = {I_dl_magnitude}·(cos({theta_s_deg}°)) = {I_dl_y:.4f} A·m")
    print(f"  I·dl_z = 0 A·m")

    print("\nStep 4: z-component of the cross product (I·dl × R)_z")
    print(f"  (I·dl × R)_z = I·dl_x·R_y - I·dl_y·R_x")
    print(f"  (I·dl × R)_z = ({I_dl_x:.4f})·({R_y:.4f}) - ({I_dl_y:.4f})·({R_x:.4f}) = {cross_product_z:.4f} A·m²")

    print("\nStep 5: Final Calculation of B_z")
    print(f"  B_z = (μ₀ / 4π) * (I·dl × R)_z / |R|³")
    print(f"  B_z = ({mu_0_div_4pi:.1e}) * ({cross_product_z:.4f}) / ({R_magnitude:.4f})³")
    print(f"  B_z = ({mu_0_div_4pi:.1e}) * ({cross_product_z:.4f}) / ({R_magnitude**3:.4f})")
    print(f"  B_z = {B_z:.4e} T")

    print(f"\nThe magnitude of the longitudinal component of the magnetic field B_z is {abs(B_z):.4e} T.")

if __name__ == '__main__':
    calculate_magnetic_field_component()