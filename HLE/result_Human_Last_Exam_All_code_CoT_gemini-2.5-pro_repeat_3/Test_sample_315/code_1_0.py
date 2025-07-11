import math

def analyze_optical_rotation():
    """
    This function calculates the rotation of polarized light for different colors
    in a D-glucose solution and explains the resulting visual phenomenon.
    """

    # --- Setup Parameters ---
    # Path length in meters, converted to decimeters (dm) for the formula
    path_length_m = 1.0
    path_length_dm = path_length_m * 10

    # Concentration of D-glucose solution in g/mL (e.g., 250g in 1L of water)
    concentration_g_ml = 0.25

    # Specific rotation for D-glucose at the sodium D-line (yellow light)
    # Unit: degrees * dm^-1 * (g/mL)^-1
    specific_rotation_yellow = 52.7

    # Wavelengths for different colors in nanometers (nm)
    lambda_red_nm = 650
    lambda_yellow_nm = 589
    lambda_blue_nm = 450

    print("--- Analysis of Optical Rotation in a Glucose Tube ---\n")
    print("Step 1: Calculate the total rotation for yellow light.")
    # Formula: total_rotation = specific_rotation * path_length * concentration
    total_rotation_yellow = specific_rotation_yellow * path_length_dm * concentration_g_ml
    print(f"The total rotation for yellow light (wavelength = {lambda_yellow_nm} nm) is:")
    print(f"α_yellow = [α] * l * c = {specific_rotation_yellow} * {path_length_dm} * {concentration_g_ml} = {total_rotation_yellow:.2f} degrees")
    print(f"This is equivalent to {total_rotation_yellow/360:.2f} full rotations.\n")


    print("Step 2: Estimate the rotation for red and blue light.")
    print("Optical rotation is approximately proportional to 1/λ², where λ is the wavelength.")
    # Estimate rotation for red light
    total_rotation_red = total_rotation_yellow * (lambda_yellow_nm / lambda_red_nm)**2
    print(f"The estimated total rotation for red light (wavelength = {lambda_red_nm} nm) is:")
    print(f"α_red = α_yellow * (λ_yellow / λ_red)² = {total_rotation_yellow:.2f} * ({lambda_yellow_nm} / {lambda_red_nm})² = {total_rotation_red:.2f} degrees")
    print(f"This is equivalent to {total_rotation_red/360:.2f} full rotations.\n")

    # Estimate rotation for blue light
    total_rotation_blue = total_rotation_yellow * (lambda_yellow_nm / lambda_blue_nm)**2
    print(f"The estimated total rotation for blue light (wavelength = {lambda_blue_nm} nm) is:")
    print(f"α_blue = α_yellow * (λ_yellow / λ_blue)² = {total_rotation_yellow:.2f} * ({lambda_yellow_nm} / {lambda_blue_nm})² = {total_rotation_blue:.2f} degrees")
    print(f"This is equivalent to {total_rotation_blue/360:.2f} full rotations.\n")

    print("--- Conclusion ---")
    print("The calculations show that blue light is rotated significantly more than red light.")
    print("As polarized white light travels down the tube, the planes of polarization for each color separate.")
    print("An observer viewing from the side sees scattered light. The color seen depends on the polarization angle at that point.")
    print("This continuous rotation, which varies by color, creates a rainbow gradient that progresses in a spiral down the tube.")

# Run the analysis
analyze_optical_rotation()