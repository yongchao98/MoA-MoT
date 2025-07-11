import math

def calculate_and_explain_sugar_tube_effect():
    """
    This script demonstrates the physics behind the optical activity of a sugar solution,
    which causes a spiraling rainbow effect.
    """

    # --- Problem Parameters & Physical Constants ---
    # These values are based on the problem description and known physical properties.
    # The 6cm diameter is not needed for this calculation.
    # The polarized filter at the end is for viewing the effect head-on;
    # we are concerned with the side view, which is due to scattering.

    tube_length_m = 1.0
    # Convert length to decimeters (dm) as it's the standard unit for the specific rotation formula.
    # path_length (L) = 1 m = 10 dm
    path_length_dm = tube_length_m * 10

    # We will assume a reasonable concentration for the D-glucose solution, e.g., 250 grams per liter.
    # concentration (c) in the formula is in g/mL.
    # 250 g/L = 0.25 g/mL
    concentration_g_per_mL = 0.25

    print("This simulation calculates the rotation angle of polarized light for different colors")
    print("at the end of a 1-meter tube of D-glucose solution.")
    print("The fact that each color rotates by a different amount is the reason for the rainbow effect.")
    print("-" * 60)
    print("Calculation Inputs:")
    print(f"Path Length (L): {path_length_dm} dm")
    print(f"Concentration (c): {concentration_g_per_mL} g/mL")

    # --- Optical Rotatory Dispersion (ORD) ---
    # Specific rotation, [α], depends on wavelength (λ). A common approximation is Biot's Law: [α] ~ 1/λ^2.
    # We can find the constant of proportionality using the known specific rotation of D-glucose
    # for the sodium D-line (yellow light).
    # [α] for λ=589.3 nm is +52.7 degrees·mL·g⁻¹·dm⁻¹.

    # Reference values
    specific_rotation_ref = 52.7  # degrees·mL·g⁻¹·dm⁻¹
    lambda_ref_nm = 589.3  # nm

    # Calculate proportionality constant k where [α] = k / λ^2
    k_constant = specific_rotation_ref * (lambda_ref_nm ** 2)

    # Define wavelengths for different representative colors (in nm)
    colors = {
        "Red": 650,
        "Yellow": 589.3,
        "Green": 540,
        "Blue": 470,
        "Violet": 420
    }
    
    print("\nFinal Rotation Angles (α = [α] * L * c):")
    print("-" * 60)
    # The final equation is α = (k / λ²) * L * c
    # Let's output each number for one example (Red light)
    red_wavelength = colors["Red"]
    red_specific_rotation = k_constant / (red_wavelength ** 2)
    red_final_angle = red_specific_rotation * path_length_dm * concentration_g_per_mL
    print(f"Example for Red light (λ={red_wavelength} nm):")
    print(f"α = ({red_specific_rotation:.2f} deg·mL·g⁻¹·dm⁻¹) * ({path_length_dm} dm) * ({concentration_g_per_mL} g/mL) = {red_final_angle:.2f} degrees")
    print("-" * 60)
    
    print(f"{'Color':<10} | {'Wavelength (nm)':<18} | {'Final Rotation (deg)':<25}")
    print("-" * 60)

    for color, wavelength in colors.items():
        # 1. Calculate specific rotation [α] for this wavelength
        specific_rotation = k_constant / (wavelength ** 2)
        # 2. Calculate the final rotation angle α at the end of the tube
        final_angle = specific_rotation * path_length_dm * concentration_g_per_mL

        print(f"{color:<10} | {wavelength:<18} | {final_angle:<25.2f}")

    print("\n--- Conclusion ---")
    print("As shown by the calculations, after traveling 1 meter, the polarization plane for Violet light")
    print("has rotated significantly more than for Red light. At every point along the tube,")
    print("the different colors have different polarization orientations. This causes a changing")
    print("color pattern when viewed from the side, which traces a spiral down the tube.")

calculate_and_explain_sugar_tube_effect()