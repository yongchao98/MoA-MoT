import numpy as np

def solve_pandora_case():
    """
    Analyzes the measurements of the star Pandora based on Planck's Law.
    """
    # Define physical constants
    h = 6.62607015e-34  # Planck constant in J·s
    c = 299792458       # Speed of light in m/s
    k = 1.380649e-23    # Boltzmann constant in J/K

    # The user's measurements
    B_measured = 1.2e15  # Spectral radiance in W·m⁻²·sr⁻¹·m⁻¹
    L_measured_nm = 400.0 # Wavelength in nm
    T_measured = 9000.0   # Temperature in K

    # Convert wavelength to SI units (meters) for the calculation
    L_measured_m = L_measured_nm * 1e-9

    print("Analyzing measurements using Planck's Law for black-body radiation.")
    print("The key insight is that the measured spectral radiance (B) is too high for the measured temperature (T).")
    print("A body at 9000 K cannot physically produce such high radiance at any wavelength.")
    print("Therefore, we conclude the temperature measurement is the most likely source of error.\n")
    print("We will now calculate the correct temperature assuming B and L are correct.")

    # --- Calculation of corrected Temperature ---
    print("The formula to solve for temperature (T) is derived from Planck's law:")
    print("T = (h * c) / (L * k * ln((2 * h * c**2) / (B * L**5) + 1))\n")

    print("The values used in the equation are:")
    print(f"B (spectral radiance) = {B_measured:.1e} W·m⁻²·sr⁻¹·m⁻¹")
    print(f"L (wavelength)        = {L_measured_m:.1e} m ({int(L_measured_nm)} nm)")
    print(f"h (Planck constant)   = {h:.4e} J·s")
    print(f"c (speed of light)    = {c:.4e} m/s")
    print(f"k (Boltzmann const)   = {k:.4e} J/K\n")

    # Perform the calculation for T
    try:
        term_in_log = (2 * h * c**2) / (B_measured * L_measured_m**5) + 1
        if term_in_log <= 1:
            print("Error: The values result in a non-positive logarithm, which is physically impossible.")
            return

        T_calculated = (h * c) / (L_measured_m * k * np.log(term_in_log))

        print("The final equation with all numbers substituted is:")
        # For readability, intermediate terms are calculated here.
        hc = h*c
        two_h_c_sq = 2 * h * c**2
        L_pow_5 = L_measured_m**5
        B_L_pow_5 = B_measured * L_pow_5

        print(f"T = ({hc:.4e}) / ({L_measured_m:.1e} * {k:.4e} * ln(({two_h_c_sq:.4e}) / ({B_L_pow_5:.4e}) + 1))")
        
        print(f"\nThe calculated correct temperature is: {T_calculated:.0f} K")

    except Exception as e:
        print(f"An error occurred during calculation: {e}")

solve_pandora_case()