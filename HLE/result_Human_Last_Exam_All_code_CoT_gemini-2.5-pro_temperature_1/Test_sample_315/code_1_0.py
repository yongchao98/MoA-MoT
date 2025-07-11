import math

def calculate_optical_rotation_distance():
    """
    This script calculates the distance along a D-glucose tube
    at which the plane of polarization for different colors of light
    is rotated by a target angle, illustrating optical rotatory dispersion.
    """
    # --- Constants for D-glucose solution ---

    # Concentration of a D-glucose solution (g/mL).
    # A moderately concentrated solution is used for a strong effect.
    c = 0.5  # g/mL

    # Constants for Drude's one-term equation for D-glucose, which models
    # how specific rotation changes with wavelength.
    # [α] = A / (λ^2 - λ_0^2)
    A = 12.33e6  # deg·nm^2·mL·g^-1·dm^-1
    lambda_0 = 150.0  # nm

    # The target angle of rotation we are calculating for.
    target_rotation = 90.0  # degrees

    # Wavelengths for different representative colors (in nanometers).
    colors = {
        "Red": 650,
        "Orange": 600,
        "Yellow": 589,
        "Green": 530,
        "Blue": 470,
        "Violet": 420
    }

    print("--- Calculation of Distance for 90° Polarization Rotation ---")
    print(f"Using a D-glucose solution with concentration c = {c} g/mL.\n")
    print("The distance is calculated using the formula: Distance(cm) = (α / ([α] * c)) * 10")
    print("where α is the target rotation, and [α] is the specific rotation for that color.\n")

    for color, wavelength in colors.items():
        # Calculate specific rotation [α] for this wavelength using Drude's equation.
        specific_rotation = A / (wavelength**2 - lambda_0**2)

        # Calculate the required path length 'l' in centimeters.
        # The base formula gives path length in decimeters (dm), so we multiply by 10.
        path_length_cm = (target_rotation / (specific_rotation * c)) * 10

        print(f"For {color} light ({wavelength} nm):")
        print(f"  Specific Rotation [α] = {specific_rotation:.2f} deg·mL·g⁻¹·dm⁻¹")
        # Final equation with all numbers shown
        print(f"  Distance (cm) = ({target_rotation} / ({specific_rotation:.2f} * {c})) * 10 = {path_length_cm:.1f} cm\n")

    print("The results show that shorter wavelengths (Violet, Blue) are rotated 90° over a")
    print("shorter distance than longer wavelengths (Red). This separation of colors along")
    print("the tube, combined with the rotational nature of the effect, creates the")
    print("appearance of a spiraling rainbow.")

calculate_optical_rotation_distance()