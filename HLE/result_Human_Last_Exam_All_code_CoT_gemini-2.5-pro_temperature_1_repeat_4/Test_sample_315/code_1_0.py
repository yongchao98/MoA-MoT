import math

def calculate_and_explain_optical_rotation():
    """
    This function calculates and explains the phenomenon of optical rotation in a glucose solution,
    leading to the appearance of a spiral rainbow.
    """
    # Plan:
    # 1. Define physical parameters for the experiment.
    # 2. Use known optical properties of D-glucose (specific rotation is wavelength-dependent).
    # 3. Calculate the total rotation and the "pitch" of the spiral for different colors.
    # 4. Print the results to demonstrate why a spiral rainbow forms.

    # Step 1: Define physical parameters
    # The problem states a 1m tall tube. In polarimetry, path length (l) is in decimeters (dm).
    # 1 m = 10 dm
    path_length_dm = 10.0

    # The concentration (c) is not given, so we'll assume a high concentration
    # typical for this demonstration, in grams per milliliter (g/mL).
    concentration_g_per_mL = 0.5

    # Step 2: Define specific rotation [α] for D-glucose for different colors.
    # This effect is called optical rotatory dispersion. Shorter wavelengths (blue) rotate more.
    # The units are degrees * mL / (g * dm). These are representative values.
    alpha_red = 45.0
    alpha_green = 60.0
    alpha_blue = 90.0

    print("Analyzing the optical rotation of D-glucose for different colors of light.")
    print(f"Parameters: Path Length = {path_length_dm} dm, Concentration = {concentration_g_per_mL} g/mL\n")

    # Step 3 & 4: Calculate total rotation and the spiral pitch for each color.
    # The pitch is the distance along the tube required for a color's polarization plane to rotate a full 360 degrees.
    # Formula for pitch (m): pitch = 360 / ([α] * 10 * c)
    
    print("--- Calculation of Spiral Pitch for Each Color ---")
    print("The 'pitch' of the spiral is the distance for one full 360-degree rotation.")
    print("A different pitch for each color results in a spiral rainbow.\n")
    
    # Red Light
    pitch_red_m = 360 / (alpha_red * 10 * concentration_g_per_mL)
    print("Red Light Pitch Calculation:")
    print(f"    360 / ({alpha_red} * 10 * {concentration_g_per_mL}) = {pitch_red_m:.2f} meters")

    # Green Light
    pitch_green_m = 360 / (alpha_green * 10 * concentration_g_per_mL)
    print("Green Light Pitch Calculation:")
    print(f"    360 / ({alpha_green} * 10 * {concentration_g_per_mL}) = {pitch_green_m:.2f} meters")

    # Blue Light
    pitch_blue_m = 360 / (alpha_blue * 10 * concentration_g_per_mL)
    print("Blue Light Pitch Calculation:")
    print(f"    360 / ({alpha_blue} * 10 * {concentration_g_per_mL}) = {pitch_blue_m:.2f} meters")

    print("\nConclusion: Because each color has a different spiral pitch, they trace separate intertwined helices down the tube, creating the visual effect of a rainbow in a spiral.")

calculate_and_explain_optical_rotation()