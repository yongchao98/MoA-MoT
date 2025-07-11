def calculate_optical_rotation():
    """
    Calculates and displays the total optical rotation for different wavelengths of light
    in a D-glucose solution, explaining the resulting visual phenomenon.
    """
    # Parameters of the experiment
    # Path length in meters, converted to decimeters (1 m = 10 dm)
    path_length_m = 1.0
    path_length_dm = path_length_m * 10

    # Assumed concentration of D-glucose solution in g/mL for a strong effect.
    concentration_g_mL = 0.5

    # Specific rotation [α] for D-glucose is dependent on the wavelength of light.
    # Shorter wavelengths are rotated more. Values are in degrees·mL·g⁻¹·dm⁻¹.
    # Data is for (Color, Wavelength in nm, Specific Rotation)
    optical_data = {
        "Red": {"wavelength": 656, "specific_rotation": 42.0},
        "Yellow (Sodium D-line)": {"wavelength": 589, "specific_rotation": 52.7},
        "Green (Mercury line)": {"wavelength": 546, "specific_rotation": 62.0},
        "Blue (Mercury line)": {"wavelength": 436, "specific_rotation": 102.0}
    }

    print("The total optical rotation (α) is calculated by the formula:")
    print("α = [α] * l * c")
    print(f"Where path length (l) = {path_length_dm} dm and concentration (c) = {concentration_g_mL} g/mL.\n")
    print("-" * 60)

    for color, data in optical_data.items():
        specific_rotation = data["specific_rotation"]
        # Calculate the total observed rotation
        total_rotation = specific_rotation * path_length_dm * concentration_g_mL
        
        print(f"For {color} light:")
        # Final equation output with all numbers
        print(f"Total Rotation α = {specific_rotation}°·mL·g⁻¹·dm⁻¹ * {path_length_dm} dm * {concentration_g_mL} g/mL = {total_rotation:.1f}°")
        print("-" * 60)
        
    print("\nAs shown by the calculations, blue light is rotated significantly more than red light.")
    print("This difference in rotation for each color causes the planes of polarization to fan out,")
    print("resulting in a visible spiral rainbow pattern along the tube.")

calculate_optical_rotation()
