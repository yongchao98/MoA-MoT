import math

def calculate_and_print_rotation():
    """
    Calculates and prints the rotation of polarization for different colors of light
    in a D-glucose solution, demonstrating the principle of optical rotatory dispersion.
    """
    # Physical parameters
    tube_length_m = 1.0  # Length of the tube in meters
    concentration_g_mL = 0.5  # Concentration of D-glucose in g/mL

    # Wavelengths for Red, Green, Blue light in nanometers (nm)
    lambda_red_nm = 650
    lambda_green_nm = 550
    lambda_blue_nm = 450

    # We can model specific rotation [α] using an approximate formula based on
    # optical rotatory dispersion: [α] ≈ k / (λ^2).
    # We find 'k' using the known value for D-glucose at the sodium D-line:
    # [α] is +52.7 deg·dm⁻¹·(g/mL)⁻¹ at λ = 589 nm.
    # k = 52.7 * (589^2) ≈ 18,345,745 deg·nm²·dm⁻¹·(g/mL)⁻¹
    k_constant = 18345745

    # Calculate specific rotation [α] for each color in deg·dm⁻¹·(g/mL)⁻¹
    spec_rot_red = k_constant / (lambda_red_nm ** 2)
    spec_rot_green = k_constant / (lambda_green_nm ** 2)
    spec_rot_blue = k_constant / (lambda_blue_nm ** 2)

    print("This script demonstrates why a rainbow spiral appears in the glucose tube.")
    print("It calculates the different rotation angles for different colors of light.\n")
    print(f"Physical setup:")
    print(f"Tube Length: {tube_length_m} m")
    print(f"Glucose Concentration: {concentration_g_mL} g/mL\n")

    print("Calculated specific rotation [α] for each color:")
    # The final equation is: Angle = [α] * length * concentration
    print(f"Red   ({lambda_red_nm} nm): [α] = {k_constant} / {lambda_red_nm}^2 = {spec_rot_red:.2f} deg·dm⁻¹·(g/mL)⁻¹")
    print(f"Green ({lambda_green_nm} nm): [α] = {k_constant} / {lambda_green_nm}^2 = {spec_rot_green:.2f} deg·dm⁻¹·(g/mL)⁻¹")
    print(f"Blue  ({lambda_blue_nm} nm): [α] = {k_constant} / {lambda_blue_nm}^2 = {spec_rot_blue:.2f} deg·dm⁻¹·(g/mL)⁻¹\n")

    print("--- Total Rotation Angle at Full Tube Length (1.0 m) ---")
    
    # Path length 'l' in decimeters (dm)
    l_dm = tube_length_m * 10
    
    # Calculate total rotation angle: α = [α] * l * c
    angle_red = spec_rot_red * l_dm * concentration_g_mL
    angle_green = spec_rot_green * l_dm * concentration_g_mL
    angle_blue = spec_rot_blue * l_dm * concentration_g_mL

    # Print the final equation with all numbers for each color
    print("Final Equation: Angle = [Specific Rotation] * [Path Length (dm)] * [Concentration (g/mL)]\n")
    print(f"Red:   Angle = {spec_rot_red:.2f} * {l_dm:.1f} * {concentration_g_mL:.1f} = {angle_red:.1f} degrees")
    print(f"Green: Angle = {spec_rot_green:.2f} * {l_dm:.1f} * {concentration_g_mL:.1f} = {angle_green:.1f} degrees")
    print(f"Blue:  Angle = {spec_rot_blue:.2f} * {l_dm:.1f} * {concentration_g_mL:.1f} = {angle_blue:.1f} degrees")

    print("\nAs you can see, over the 1m length, the colors have rotated by very different amounts.")
    print("This continuous, color-dependent rotation creates the appearance of a spiral rainbow.")

calculate_and_print_rotation()