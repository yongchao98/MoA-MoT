import math

def calculate_rotation():
    """
    This function calculates and demonstrates the principle of optical rotatory dispersion
    for a D-glucose solution, which explains the observed visual phenomenon.
    """
    # The formula for optical rotation is:
    # θ = [α] * L * c
    # where:
    # θ is the angle of rotation in degrees.
    # [α] is the specific rotation, which depends on the wavelength (color) of light.
    # L is the path length in decimeters (dm).
    # c is the concentration in g/mL.

    # --- Setup Parameters ---
    # The tube is 1m tall.
    path_length_m = 1.0
    # Convert path length to decimeters for the formula.
    path_length_dm = path_length_m * 10

    # We assume a reasonable concentration for a strong effect.
    # Let's say 200 grams of glucose per liter of water.
    concentration_g_L = 200.0
    # Convert concentration to g/mL for the formula.
    concentration_g_ml = concentration_g_L / 1000.0

    # Specific rotation [α] for D-glucose varies with wavelength (color).
    # We use known values for different spectral lines as representatives for colors.
    # Units are in degrees*mL/(g*dm).
    specific_rotation = {
        'Violet (436 nm)': 102.0,
        'Green (546 nm)': 62.0,
        'Yellow (589 nm)': 52.7
    }

    print("The total rotation of the plane of polarization for different colors of light")
    print(f"after traveling {path_length_m} m through a D-glucose solution is calculated below.")
    print("-" * 70)

    for color, alpha in specific_rotation.items():
        # Calculate the total rotation angle
        total_rotation = alpha * path_length_dm * concentration_g_ml

        # Print the equation with the numbers substituted in
        print(f"For {color}:")
        print(f"  Rotation Angle θ = [α] * L * c")
        print(f"  θ = {alpha} deg*mL/(g*dm) * {path_length_dm} dm * {concentration_g_ml} g/mL")
        print(f"  θ = {total_rotation:.1f} degrees\n")

    print("-" * 70)
    print("This difference in rotation angles for each color is called optical rotatory dispersion.")
    print("It causes the polarization planes to separate, leading to a spiral rainbow of scattered light visible from the side.")

calculate_rotation()