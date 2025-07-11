def calculate_optical_rotation():
    """
    Calculates the total optical rotation for different wavelengths of light
    passing through a D-glucose solution, demonstrating the principle of
    optical rotatory dispersion.
    """
    # Parameters of the experiment
    path_length_m = 1.0  # Path length in meters
    path_length_dm = path_length_m * 10  # Convert meters to decimeters (dm) for the formula
    concentration_g_cm3 = 0.25  # Assumed concentration in g/cm^3 (e.g., 25g per 100mL)

    # Specific rotation [α] for D-glucose varies with wavelength (λ).
    # Units are in degrees * cm^3 / (g * dm).
    # These are approximate literature values.
    specific_rotations = {
        'Red (650 nm)': 46.0,
        'Yellow (589 nm)': 52.7,
        'Green (546 nm)': 62.3,
        'Blue (450 nm)': 95.0,
        'Violet (436 nm)': 101.4,
    }

    print(f"Calculating total optical rotation for a {path_length_m}m tube of D-glucose solution.\n")
    print("The formula is: Observed Rotation = Specific Rotation * Concentration * Path Length\n")

    for color, spec_rot in specific_rotations.items():
        # Calculate the observed rotation using the formula: α = [α] * l * c
        observed_rotation = spec_rot * path_length_dm * concentration_g_cm3

        print(f"For {color}:")
        # The prompt requires outputting each number in the final equation.
        print(f"  Observed Rotation = {spec_rot} deg*cm^3/(g*dm) * {concentration_g_cm3} g/cm^3 * {path_length_dm} dm")
        print(f"  Result: {observed_rotation:.2f} degrees of rotation.\n")

    print("As you can see, shorter wavelengths (like blue/violet) rotate significantly more than longer")
    print("wavelengths (like red). This separation of rotation angles causes the spiraling rainbow effect.")

calculate_optical_rotation()
