import math

def calculate_optical_rotation():
    """
    This script demonstrates optical rotatory dispersion for a D-glucose solution.
    It calculates the total rotation angle for different colors of light passing
    through a 1-meter tube, which explains the resulting spiral rainbow effect.
    """
    # Plan:
    # 1. Define the physical parameters of the experiment based on the problem description.
    # 2. Use representative data for the specific rotation of D-glucose at different wavelengths (colors).
    # 3. Apply Biot's Law to calculate the total rotation for each color.
    # 4. Print the calculation for each color, showing how shorter wavelengths are rotated more.

    # Step 1: Define physical parameters
    # The tube is 1m tall, so the path length 'l' is 1.0 m.
    # Biot's law uses decimeters (dm), so l = 10 dm.
    path_length_dm = 10.0

    # The concentration 'c' is not specified. We'll assume a standard
    # concentration for demonstrations, e.g., 200 grams per Liter, which is 0.2 g/mL.
    concentration_g_ml = 0.2

    print(f"Analyzing optical rotation in a {path_length_dm} dm tube of D-glucose solution with a concentration of {concentration_g_ml} g/mL.")
    print("The final appearance from the side is due to optical rotatory dispersion, where rotation angle depends on wavelength.")
    print("-" * 80)

    # Step 2: Define specific rotation data for different wavelengths.
    # The specific rotation [a] is a property of the substance. Units are deg*mL/(g*dm).
    # These are typical values for D-glucose, showing that [a] increases as wavelength decreases.
    color_data = {
        "Red (656 nm)": 46.8,
        "Yellow (589 nm)": 52.7,
        "Green (546 nm)": 63.5,
        "Blue (436 nm)": 101.4,
    }

    # Step 3 & 4: Calculate and print the total rotation for each color.
    print("Formula: Total Rotation = [Specific Rotation] * Path Length * Concentration\n")

    for color, specific_rotation in color_data.items():
        # Calculate total rotation using Biot's Law
        total_rotation = specific_rotation * path_length_dm * concentration_g_ml

        # Print the equation with the values substituted in
        print(f"For {color}:")
        print(f"    Total Rotation = {specific_rotation} * {path_length_dm} * {concentration_g_ml}")
        print(f"    Result = {total_rotation:.1f} degrees of rotation\n")

    print("-" * 80)
    print("Since each color's plane of polarization rotates by a different amount, the color")
    print("scattered towards a side observer changes along the tube's length, creating a spiral rainbow.")

calculate_optical_rotation()

<<<C>>>