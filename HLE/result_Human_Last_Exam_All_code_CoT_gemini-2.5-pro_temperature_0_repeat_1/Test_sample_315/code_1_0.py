import math

def solve_optical_rotation():
    """
    This script calculates the optical rotation for different colors of light
    in a D-glucose solution to explain the visual phenomenon described.
    """
    # Plan:
    # 1. Define the physical parameters of the experiment based on the problem description.
    # 2. Use known approximate values for the specific rotation of D-glucose for different
    #    wavelengths (colors) to model Optical Rotatory Dispersion.
    # 3. Apply the formula for optical rotation to find the total rotation angle for each color.
    # 4. Print the results and explain how they lead to the final answer.

    # Step 1: Define parameters
    # Path length (l) is 1m. The standard unit for the formula is decimeters (dm).
    # 1 m = 10 dm
    path_length_dm = 10.0

    # Let's assume a reasonable concentration (c) for the solution, e.g., 250 g/L.
    # The standard unit is g/mL. 250 g/L = 0.25 g/mL
    concentration_g_ml = 0.25

    # Step 2: Define specific rotation [α] for D-glucose in (degrees * mL) / (g * dm).
    # Specific rotation depends on wavelength (color). Shorter wavelengths are rotated more.
    # These are representative values.
    specific_rotation_red = 40.0   # For red light (~650 nm)
    specific_rotation_green = 60.0 # For green light (~550 nm)
    specific_rotation_blue = 110.0  # For blue light (~450 nm)

    # Step 3: Calculate the total rotation angle (theta) for each color.
    # Formula: theta = specific_rotation * concentration * path_length
    theta_red = specific_rotation_red * concentration_g_ml * path_length_dm
    theta_green = specific_rotation_green * concentration_g_ml * path_length_dm
    theta_blue = specific_rotation_blue * concentration_g_ml * path_length_dm

    # Step 4: Print the explanation and results.
    print("--- Optical Rotation Calculation for D-Glucose Solution ---")
    print("The appearance of the tube is due to optical activity combined with optical rotatory dispersion.")
    print("This means the plane of polarized light is rotated by the glucose, and the angle of rotation depends on the light's color.")
    print("\nLet's calculate the total rotation for different colors over the 1m path length:")
    print(f"\nParameters:")
    print(f"  Path Length: {path_length_dm} dm")
    print(f"  Concentration: {concentration_g_ml} g/mL")

    print("\nCalculation: Angle = [α] * c * l")
    print("-" * 50)
    print(f"Rotation for RED light:   {specific_rotation_red:.1f} * {concentration_g_ml} * {path_length_dm} = {theta_red:.1f} degrees")
    print(f"Rotation for GREEN light: {specific_rotation_green:.1f} * {concentration_g_ml} * {path_length_dm} = {theta_green:.1f} degrees")
    print(f"Rotation for BLUE light:  {specific_rotation_blue:.1f} * {concentration_g_ml} * {path_length_dm} = {theta_blue:.1f} degrees")
    print("-" * 50)

    print("\nConclusion:")
    print("The final polarization angles for red, green, and blue light are separated by a large amount.")
    print("This rotation happens continuously along the tube. The combination of continuous rotation and color separation creates a rainbow.")
    print("Because the rotation happens around the axis of the tube, the rainbow pattern appears as a spiral.")

solve_optical_rotation()