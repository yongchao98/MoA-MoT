import math

def calculate_optical_rotation():
    """
    This script demonstrates the principle of optical rotatory dispersion in a D-glucose solution.
    It calculates the rotation of the plane of polarized light for different colors (wavelengths)
    at various distances along a tube. This illustrates why a spiral rainbow effect is observed.
    """

    # --- Parameters ---
    # Length of the tube in meters.
    tube_length_m = 1.0
    # Concentration of D-glucose solution in grams per milliliter.
    # A 25% w/v solution is reasonably concentrated for a strong effect.
    concentration_g_per_ml = 0.25

    # Specific rotation [alpha] in (deg * mL) / (g * dm) for different wavelengths of light.
    # The value is dependent on the wavelength (color), which is the principle of optical rotatory dispersion.
    # These are typical approximate values for D-glucose.
    specific_rotations = {
        'Red (~650nm)': 40.0,
        'Green (~546nm)': 63.0,
        'Blue (~436nm)': 100.0
    }

    # --- Calculation and Output ---
    print("Demonstration of Optical Rotation in D-Glucose Solution\n")
    print(f"Tube Length: {tube_length_m} m")
    print(f"Concentration (c): {concentration_g_per_ml} g/mL")
    print("\nThe angle of rotation (θ) is calculated using Biot's Law: θ = [α] * l * c")
    print("where [α] is the specific rotation, l is the path length, and c is the concentration.")
    print("-" * 75)

    # We will check the rotation at several points along the tube.
    distances_to_check_m = [0.0, 0.25, 0.5, 0.75, 1.0]

    for distance_m in distances_to_check_m:
        # Biot's law requires path length (l) in decimeters (dm). 1 m = 10 dm.
        path_length_dm = distance_m * 10
        print(f"\nAt a distance of {distance_m} m (path length l = {path_length_dm:.1f} dm):")

        for color, alpha in specific_rotations.items():
            # This is the core calculation for optical rotation.
            rotation_angle = alpha * path_length_dm * concentration_g_per_ml
            
            # We output the full equation with the numbers plugged in.
            equation_str = f"θ = {alpha:.1f} * {path_length_dm:.1f} * {concentration_g_per_ml}"
            result_str = f"{rotation_angle:.1f} degrees"
            print(f"  - For {color.ljust(15)}: {equation_str.ljust(35)} = {result_str}")

    print("\n" + "-" * 75)
    print("Conclusion: At any point > 0m, the rotation angles are different for each color.")
    print("This separation of colors by rotation angle creates a rainbow.")
    print("Because the rotation happens around the axis of light travel, the pattern appears as a spiral.")

calculate_optical_rotation()