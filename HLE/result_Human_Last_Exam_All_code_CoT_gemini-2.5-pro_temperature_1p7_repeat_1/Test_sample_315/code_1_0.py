import math

def calculate_optical_rotation():
    """
    This function demonstrates the principle of optical rotatory dispersion
    for a D-glucose solution, which explains the observed phenomenon.
    """
    
    # --- Parameters of the setup ---
    # Path length in meters
    length_m = 1.0
    # The formula for specific rotation uses path length in decimeters (dm)
    length_dm = length_m * 10
    
    # We assume a reasonable concentration for a strong effect.
    # The standard unit is g/mL. (e.g., 25 g in 100 mL is 0.25 g/mL)
    concentration_g_mL = 0.50

    print(f"Analysis for a {length_m}m tube with D-glucose solution at {concentration_g_mL} g/mL.\n")

    # Specific rotation [α] for D-glucose is dependent on the wavelength (color) of light.
    # Data is for the D-line of sodium, and green and violet lines of mercury.
    # The format is { 'Color (Wavelength nm)': specific_rotation_angle }
    specific_rotations = {
        'Violet (436 nm)': 102.7,
        'Green (546 nm)': 62.0,
        'Yellow (589 nm)': 52.7,
    }

    print("Equation: Total Rotation = [α] * L * c")
    print(f"Where L = {length_dm} dm and c = {concentration_g_mL} g/mL\n")

    # Calculate and print the total rotation for each color
    for color, spec_rot in specific_rotations.items():
        total_rotation = spec_rot * length_dm * concentration_g_mL
        
        # We present the final calculation as the required equation
        print(f"Total Rotation for {color}:")
        print(f"{total_rotation:.1f} degrees = {spec_rot} * {length_dm} * {concentration_g_mL}")
        print(f"(This is equivalent to {total_rotation/360:.2f} full turns)\n")

    print("Since each color rotates by a different amount, the color of light scattered")
    print("to a side-viewer changes along the tube's length, creating a spiral rainbow.")

calculate_optical_rotation()
