import math

def calculate_and_explain_optical_rotation():
    """
    Calculates and explains the optical rotation of different colors of light
    in a D-glucose solution, which leads to a spiral rainbow effect.
    """

    # --- Problem Parameters ---
    # The path length (l) of the tube is 1 meter. For Biot's Law, we convert this
    # to decimeters (dm), so 1 m = 10 dm.
    path_length_dm = 10.0

    # Let's assume a reasonable concentration (c) for the solution, for example,
    # 250 g/L. This is converted to g/mL for the formula, so 250 g/L = 0.25 g/mL.
    concentration_g_per_ml = 0.25

    # Specific rotation [α] for D-glucose is dependent on the wavelength (color).
    # These are typical approximate values in units of degrees·mL·g⁻¹·dm⁻¹.
    # Shorter wavelengths (blue) are rotated more than longer wavelengths (red).
    specific_rotation_red = 40.0   # For red light (~650 nm)
    specific_rotation_blue = 100.0  # For blue light (~450 nm)

    print("This script demonstrates why a tube of D-glucose solution can create a spiral rainbow.")
    print("The effect is due to optical rotatory dispersion, where different colors of light rotate by different amounts.")
    print("\nWe use Biot's Law: Total Rotation (θ) = Specific Rotation ([α]) * Path Length (l) * Concentration (c)\n")

    # --- Calculations ---
    # Calculate the total rotation for red and blue light.
    total_rotation_red = specific_rotation_red * path_length_dm * concentration_g_per_ml
    total_rotation_blue = specific_rotation_blue * path_length_dm * concentration_g_per_ml

    # --- Output Results with Equation Numbers ---
    print("--- Calculation for Red Light ---")
    print(f"θ_red = [α] * l * c")
    print(f"θ_red = {specific_rotation_red} * {path_length_dm} * {concentration_g_per_ml}")
    print(f"Total Rotation for Red Light: {total_rotation_red:.1f} degrees\n")

    print("--- Calculation for Blue Light ---")
    print(f"θ_blue = [α] * l * c")
    print(f"θ_blue = {specific_rotation_blue} * {path_length_dm} * {concentration_g_per_ml}")
    print(f"Total Rotation for Blue Light: {total_rotation_blue:.1f} degrees\n")

    # --- Conclusion ---
    print("Conclusion:")
    print(f"After 1 meter, the polarization of red light rotates {total_rotation_red:.1f} degrees, while blue light rotates {total_rotation_blue:.1f} degrees.")
    print("This difference means the polarization planes for each color are constantly twisting relative to each other down the tube.")
    print("When viewed from the side, this twisting pattern of polarized colors manifests as a rainbow that spirals along the tube's length.")

if __name__ == '__main__':
    calculate_and_explain_optical_rotation()

<<<C>>>