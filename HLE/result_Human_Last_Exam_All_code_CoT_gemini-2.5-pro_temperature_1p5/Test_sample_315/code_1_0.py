import math

def calculate_and_display_rotation():
    """
    Calculates and displays the optical rotation of different colors of light
    in a D-glucose solution to explain the spiral rainbow effect.
    """
    # Plan:
    # 1. Define physical constants for the experiment.
    # 2. Model the wavelength-dependence of specific rotation [a].
    #    Specific rotation is higher for shorter wavelengths. We'll use a
    #    simple but effective model based on known data for D-glucose.
    # 3. Use Biot's Law (Angle = [a] * L * c) to calculate the rotation angle
    #    for different colors at various distances down the tube.
    # 4. Print the results to show how the polarization plane for each color
    #    rotates differently, which would appear as a spiral from the side.
    # 5. Fulfill the requirement to print each number in a final sample calculation.

    # --- Step 1: Define Constants ---
    # Tube length in decimeters (dm). Biot's Law traditionally uses dm.
    # 1 meter = 10 dm.
    L_tube_dm = 10.0
    # Concentration of D-glucose solution in g/mL. Let's assume a 20% w/v solution.
    c_glucose_g_ml = 0.20

    # --- Step 2: Model Specific Rotation [a] ---
    # The specific rotation [a] depends on wavelength.
    # Unit: degrees * mL / (g * dm)
    # A known value for D-glucose: [a] at 589 nm (Yellow) is +52.7
    # We can create a simple model: [a](lambda) = k / (lambda^2),
    # which captures the principle that shorter wavelengths are rotated more.
    # We calibrate the constant 'k' using the known value.
    specific_rotation_ref = 52.7
    lambda_ref = 589  # nm
    k_constant = specific_rotation_ref * (lambda_ref ** 2)

    def get_specific_rotation(wavelength_nm):
        """Calculates specific rotation for a given wavelength."""
        return k_constant / (wavelength_nm ** 2)

    # Wavelengths for Red, Green, and Blue light
    wavelengths = {'Red': 650, 'Green': 530, 'Blue': 470}

    print("--- Simulating Polarization Rotation in D-Glucose Solution ---")
    print(f"This simulation shows why a spiral rainbow appears.")
    print("Different colors of light rotate at different rates as they travel down the tube.")
    print("The table shows the total rotation angle for each color at different distances.\n")
    print(f"{'Distance (m)':<15} | {'Red Angle (deg)':<20} | {'Green Angle (deg)':<20} | {'Blue Angle (deg)':<20}")
    print("-" * 80)

    # --- Step 3: Calculate and Display Rotation at Different Points ---
    # Iterate through 11 steps from 0m to 1m
    for i in range(11):
        distance_m = i * 0.1
        distance_dm = distance_m * 10

        # Calculate rotation for each color at this distance
        angle_r = get_specific_rotation(wavelengths['Red']) * distance_dm * c_glucose_g_ml
        angle_g = get_specific_rotation(wavelengths['Green']) * distance_dm * c_glucose_g_ml
        angle_b = get_specific_rotation(wavelengths['Blue']) * distance_dm * c_glucose_g_ml

        print(f"{distance_m:<15.1f} | {angle_r:<20.2f} | {angle_g:<20.2f} | {angle_b:<20.2f}")

    print("-" * 80)
    print("\n--- Final Calculation Breakdown ---")
    print("Here are the numbers for the final equation for Blue light at the end of the tube (1.0m):")
    # --- Step 5: Print numbers for the final equation ---
    final_distance_dm = 10.0
    blue_wavelength = wavelengths['Blue']
    blue_spec_rot = get_specific_rotation(blue_wavelength)
    final_angle_blue = blue_spec_rot * final_distance_dm * c_glucose_g_ml
    
    print("Equation: Final Angle = Specific Rotation * Path Length * Concentration")
    print(f"Calculated Angle = {blue_spec_rot:.2f} * {final_distance_dm:.1f} * {c_glucose_g_ml:.2f} = {final_angle_blue:.2f} degrees")

# Execute the function
calculate_and_display_rotation()