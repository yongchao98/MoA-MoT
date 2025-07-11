import math

def solve_and_explain():
    """
    Calculates the rotation of light in a D-glucose solution to explain the visual appearance.
    """
    # 1. Define physical parameters from the problem and standard values.
    # The path length of the light in the tube. The formula for optical rotation
    # uses path length in decimeters (dm). 1 m = 10 dm.
    path_length_m = 1.0
    path_length_dm = path_length_m * 10

    # The concentration of the D-glucose solution is not given. We'll use a
    # plausible concentration of 0.2 g/mL for a clear demonstration.
    concentration_g_ml = 0.2

    # 2. Define specific rotation ([α]) for different wavelengths (colors).
    # Specific rotation is an intrinsic property of D-glucose and it changes with
    # the wavelength of light (Optical Rotatory Dispersion).
    # Units are in degrees / (dm * (g/mL)).
    # We use approximate values for red, yellow, and blue light.
    spec_rotation_red_approx = 45.0   # For red light (~650 nm)
    spec_rotation_yellow_Na = 52.7  # Standard value for yellow Sodium D-line (~589 nm)
    spec_rotation_blue_approx = 78.0  # For blue light (~480 nm)

    print("This problem involves optical activity, where a substance rotates the plane of polarized light.")
    print("The amount of rotation depends on the light's color (wavelength).\n")
    print(f"Physical Parameters:")
    print(f"Path Length: {path_length_m} m = {path_length_dm} dm")
    print(f"Assumed Concentration: {concentration_g_ml} g/mL\n")

    # 3. Calculate the total rotation for each color using the formula:
    # Total Rotation = [Specific Rotation] * Path Length (dm) * Concentration (g/mL)
    print("Calculating the total rotation for different colors...")
    print("Formula: Total Rotation = [α] * L * c\n")

    # Calculation for Red Light
    total_rotation_red = spec_rotation_red_approx * path_length_dm * concentration_g_ml
    print(f"Rotation for Red Light:   {spec_rotation_red_approx} * {path_length_dm} * {concentration_g_ml} = {total_rotation_red:.1f} degrees")

    # Calculation for Yellow Light
    total_rotation_yellow = spec_rotation_yellow_Na * path_length_dm * concentration_g_ml
    print(f"Rotation for Yellow Light: {spec_rotation_yellow_Na} * {path_length_dm} * {concentration_g_ml} = {total_rotation_yellow:.1f} degrees")

    # Calculation for Blue Light
    total_rotation_blue = spec_rotation_blue_approx * path_length_dm * concentration_g_ml
    print(f"Rotation for Blue Light:  {spec_rotation_blue_approx} * {path_length_dm} * {concentration_g_ml} = {total_rotation_blue:.1f} degrees\n")

    # 4. Explain the implication of the results.
    print("Conclusion:")
    print(f"After traveling 1 meter, the polarization of blue light is rotated by {total_rotation_blue:.1f} degrees,")
    print(f"while red light is rotated by only {total_rotation_red:.1f} degrees.")
    print("This difference means the polarization plane for each color is pointed in a different direction at any given point along the tube.")
    print("When viewed from the side, the scattered light we see has a color that depends on our viewing angle relative to these rotated polarization planes.")
    print("This creates a continuous, shifting rainbow of colors that appears to spiral down the tube.")

solve_and_explain()
<<<C>>>