import math

def calculate_magnetic_field_component():
    """
    Calculates the longitudinal component of the magnetic field (Bz) at a target point
    due to a current element on a cylinder surface, based on the Biot-Savart Law.
    """
    # --- Given Parameters ---
    # Permeability of free space / 4*pi
    mu0_over_4pi = 1e-7  # T·m/A

    # Target field point (Cartesian)
    x_f, y_f, z_f = 0.1, -0.1, 0.2

    # Source point (Cylindrical)
    r_s_cyl = 0.34
    theta_s_deg = 60
    z_s_cyl = 0.5

    # Azimuthal current element magnitude (interpretation of the given J_theta)
    # The problem gives a surface current density J = 40 A/m^2. To find a field in Tesla from
    # a point source, we interpret this as the magnitude of a differential current element, I*dl = 40 A·m.
    Idl_mag = 40.0  # A·m

    print("--- Problem Setup ---")
    print(f"Target field point (xf, yf, zf) = ({x_f}, {y_f}, {z_f}) m")
    print(f"Source point (r, θ, z) = ({r_s_cyl} m, {theta_s_deg}°, {z_s_cyl} m)")
    print(f"Interpreted current element magnitude |Idl| = {Idl_mag} A·m")
    print(f"Constant μ₀/(4π) = {mu0_over_4pi} T·m/A\n")

    # --- Step 1: Convert source point to Cartesian coordinates ---
    theta_s_rad = math.radians(theta_s_deg)
    x_s = r_s_cyl * math.cos(theta_s_rad)
    y_s = r_s_cyl * math.sin(theta_s_rad)
    z_s = z_s_cyl
    print("--- Step 1: Source Point in Cartesian Coordinates ---")
    print(f"Source point (xs, ys, zs) = ({x_s:.4f}, {y_s:.4f}, {z_s:.4f}) m\n")

    # --- Step 2: Calculate the displacement vector R = r_f - r_s ---
    R_x = x_f - x_s
    R_y = y_f - y_s
    R_z = z_f - z_s
    print("--- Step 2: Displacement Vector R ---")
    print(f"R = (xf-xs, yf-ys, zf-zs)")
    print(f"R = ({R_x:.4f}, {R_y:.4f}, {R_z:.4f}) m\n")
    
    # --- Step 3: Define the current element vector Idl ---
    # The direction is azimuthal (theta_hat), which in Cartesian is (-sin(theta), cos(theta), 0)
    Idl_x = -Idl_mag * math.sin(theta_s_rad)
    Idl_y = Idl_mag * math.cos(theta_s_rad)
    Idl_z = 0.0
    print("--- Step 3: Current Element Vector Idl ---")
    print(f"Idl direction is theta_hat = (-sin({theta_s_deg}°), cos({theta_s_deg}°), 0)")
    print(f"Idl = (Idl_x, Idl_y, Idl_z) = ({Idl_x:.4f}, {Idl_y:.4f}, {Idl_z:.4f}) A·m\n")

    # --- Step 4: Calculate components for Biot-Savart Law ---
    # Numerator is the z-component of the cross product (Idl x R)
    # (Idl x R)_z = Idl_x * R_y - Idl_y * R_x
    numerator = Idl_x * R_y - Idl_y * R_x
    
    # Denominator is the magnitude of R cubed, |R|^3
    R_mag_sq = R_x**2 + R_y**2 + R_z**2
    R_mag = math.sqrt(R_mag_sq)
    denominator = R_mag**3

    # --- Step 5: Final Calculation ---
    print("--- Step 4 & 5: Applying Biot-Savart Law for Bz ---")
    print("Formula: B_z = (μ₀/4π) * (Idl_x * R_y - Idl_y * R_x) / |R|³")
    print("\nPlugging in the numbers:")
    print(f"Numerator = ({Idl_x:.4f}) * ({R_y:.4f}) - ({Idl_y:.4f}) * ({R_x:.4f})")
    print(f"Numerator = {numerator:.4f} A·m²")
    print(f"Denominator |R|³ = (({R_x:.4f})² + ({R_y:.4f})² + ({R_z:.4f})²)^(3/2)")
    print(f"Denominator = {denominator:.4f} m³")
    
    # Calculate Bz
    B_z = mu0_over_4pi * numerator / denominator
    
    print("\nFinal Equation:")
    print(f"B_z = {mu0_over_4pi} * ({numerator:.4f}) / ({denominator:.4f})")
    print(f"B_z = {B_z:.4e} T")

    # The question asks for the magnitude
    magnitude_B_z = abs(B_z)
    print(f"\nThe magnitude of the longitudinal component is |B_z| = {magnitude_B_z:.4e} T")

    return magnitude_B_z

if __name__ == '__main__':
    final_answer = calculate_magnetic_field_component()
    # The final answer is expected in a specific format
    # The result is approximately 1.2017e-5 T, or 12.017 uT.
    # print(f"\n<<<{final_answer}>>>")
    # For a prettier format:
    # print(f"\n<<<{final_answer:.3e}>>>")
    # Let's provide the final result without the wrapper in the main execution block.
    # The final instruction just says return the answer.

calculate_magnetic_field_component()
# We are asked for a final answer block. Let's compute it again and format it.
B_z_val = 1.201662955561338e-05
# This is approximately 1.2017 x 10^-5
# Wrapping the final numerical value.
# I will use scientific notation for precision.
final_value_for_wrapper = B_z_val

print(f"\n<<<1.2017e-05>>>")
