import math

def calculate_and_explain_rotation():
    """
    Calculates the optical rotation for different wavelengths of light in a D-glucose solution
    to demonstrate the principle of optical rotatory dispersion.
    """
    # --- Parameters ---
    # Path length (l) in decimeters (1 m = 10 dm)
    l = 10.0
    # Assumed concentration (c) in g/mL (e.g., 200 g/L = 0.2 g/mL)
    c = 0.2
    # Wavelength for red light in nm
    lambda_red = 650.0
    # Wavelength for violet light in nm
    lambda_violet = 400.0

    # Drude's one-term equation parameters for glucose in water (approximate)
    # A is a constant in deg*nm^2*mL/(g*dm)
    A = 12.33e6
    # lambda_0 is the wavelength of the ultraviolet absorption band in nm
    lambda_0 = 150.0

    print("This script demonstrates Optical Rotatory Dispersion for D-glucose.")
    print(f"We will calculate the total rotation for red ({lambda_red} nm) and violet ({lambda_violet} nm) light.")
    print(f"Path length (l) = {l} dm, Concentration (c) = {c} g/mL.\n")

    # --- Calculation for Red Light ---
    spec_rot_red = A / (lambda_red**2 - lambda_0**2)
    obs_rot_red = spec_rot_red * l * c

    print("--- For Red Light (λ = 650 nm) ---")
    print("1. Calculate Specific Rotation [α] using Drude's equation:")
    print(f"   [α] = A / (λ^2 - λ_0^2)")
    print(f"   [α] = {A} / ({lambda_red}^2 - {lambda_0}^2)")
    print(f"   [α] = {spec_rot_red:.2f} deg*mL/(g*dm)\n")

    print("2. Calculate Observed Rotation α:")
    print(f"   α = [α] * l * c")
    print(f"   α = {spec_rot_red:.2f} * {l} * {c}")
    print(f"   Total Rotation for Red Light = {obs_rot_red:.2f} degrees\n")

    # --- Calculation for Violet Light ---
    spec_rot_violet = A / (lambda_violet**2 - lambda_0**2)
    obs_rot_violet = spec_rot_violet * l * c

    print("--- For Violet Light (λ = 400 nm) ---")
    print("1. Calculate Specific Rotation [α] using Drude's equation:")
    print(f"   [α] = A / (λ^2 - λ_0^2)")
    print(f"   [α] = {A} / ({lambda_violet}^2 - {lambda_0}^2)")
    print(f"   [α] = {spec_rot_violet:.2f} deg*mL/(g*dm)\n")

    print("2. Calculate Observed Rotation α:")
    print(f"   α = [α] * l * c")
    print(f"   α = {spec_rot_violet:.2f} * {l} * {c}")
    print(f"   Total Rotation for Violet Light = {obs_rot_violet:.2f} degrees\n")

    # --- Conclusion ---
    difference = obs_rot_violet - obs_rot_red
    print("--- Conclusion ---")
    print(f"As shown, the polarization plane for violet light rotates {obs_rot_violet:.2f} degrees,")
    print(f"while red light's plane rotates only {obs_rot_red:.2f} degrees.")
    print(f"This difference of {difference:.2f} degrees over the 1m tube is what separates the colors.")
    print("This differential rotation for all colors in white light creates the appearance of a spiral rainbow.")

if __name__ == '__main__':
    calculate_and_explain_rotation()