import math

def calculate_rotation():
    """
    Calculates the total rotation of the plane of polarization for different
    wavelengths of light passing through a D-glucose solution.
    """
    # --- Constants and Assumptions ---

    # Length of the tube in decimeters (dm). 1 m = 10 dm.
    length_dm = 10.0

    # Assumed concentration of D-glucose solution in g/mL.
    concentration_g_ml = 0.2

    # Known specific rotation for D-glucose at the sodium D-line (589 nm).
    # Units: degrees * mL / (g * dm)
    specific_rotation_589nm = 52.7
    lambda_589nm = 589.0

    # --- Wavelength-dependent rotation model ---
    # We model specific rotation [a] as [a] = k / lambda^2.
    # We find the constant k using the known value at 589 nm.
    k = specific_rotation_589nm * (lambda_589nm ** 2)

    # Wavelengths for Red, Green, and Blue light in nm.
    wavelengths = {
        'Red': 650.0,
        'Green': 550.0,
        'Blue': 450.0
    }

    print("This script calculates the total rotation angle for different colors of light")
    print("in a 1m tube of D-glucose solution based on the formula: a = [a] * l * c.")
    print(f"Using Length (l) = {length_dm} dm and Concentration (c) = {concentration_g_ml} g/mL.\n")

    # --- Calculation for each color ---
    for color, lam in wavelengths.items():
        # Calculate specific rotation for this wavelength
        specific_rotation_lambda = k / (lam ** 2)

        # Calculate total rotation angle using Biot's Law: angle = [a] * l * c
        total_rotation_angle = specific_rotation_lambda * length_dm * concentration_g_ml

        print(f"For {color} light ({int(lam)} nm):")
        # Print the equation with all the numbers
        print(f"Total Rotation = Specific Rotation * Length * Concentration")
        print(f"Total Rotation = {specific_rotation_lambda:.2f} deg*mL/(g*dm) * {length_dm} dm * {concentration_g_ml} g/mL")
        print(f"Result: The polarization plane is rotated by {total_rotation_angle:.2f} degrees.\n")

    print("The different rotation angles for each color cause the light to disperse into a spectrum,")
    print("which appears as a spiraling rainbow when viewed from the side.")


if __name__ == '__main__':
    calculate_rotation()
