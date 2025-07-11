import math

def solve_pandora_problem():
    """
    Analyzes the measurements of the white dwarf Pandora and determines if they are consistent
    with Planck's Law.
    """
    # Physical constants (using CODATA 2018 values for high precision)
    h = 6.62607015e-34  # Planck constant (J*s)
    c = 299792458      # Speed of light (m/s)
    k = 1.380649e-23     # Boltzmann constant (J/K)

    # Given measurements
    L_nm = 400.0         # Wavelength in nanometers
    L_m = L_nm * 1e-9    # Wavelength in meters
    T_k = 9000.0         # Surface temperature in Kelvin
    B_measured = 1.2e15  # Spectral radiance in W/(m^2*sr*m)

    # Step 1: Check for consistency by calculating expected spectral radiance from T and L
    # Using Planck's Law: B(L, T) = (2*h*c^2 / L^5) / (exp(h*c / (L*k*T)) - 1)
    try:
        exponent = (h * c) / (L_m * k * T_k)
        # Avoid overflow for very large exponents, though not an issue here
        if exponent > 700:
             B_expected = 0
        else:
            numerator = 2 * h * c**2
            denominator = (L_m**5) * (math.exp(exponent) - 1)
            B_expected = numerator / denominator
    except (ValueError, OverflowError):
        B_expected = float('inf')

    # Step 2: Compare measured B with expected B. A small tolerance (e.g., 5%) is used.
    # If the values were consistent, the relative difference would be small.
    relative_diff = abs(B_expected - B_measured) / B_measured
    
    if relative_diff < 0.05:
        # The result looks okay
        print("0")
    else:
        # The result is not okay. We need to find the likely error.
        # Astrophysical context: DB-class white dwarfs typically have T > 12,000 K.
        # The given T=9000K is an outlier. Thus, T is the most likely error.
        # We solve for the corrected temperature T_corr using the measured B and L.
        # From Planck's law: T = (h*c / (L*k)) / ln((2*h*c^2 / (L^5 * B)) + 1)
        
        term_inside_ln = (2 * h * c**2) / (L_m**5 * B_measured) + 1
        
        if term_inside_ln <= 1:
            print("Error: Cannot calculate temperature with the given B and L values (logarithm input is not positive).")
            return

        T_corrected = ((h * c) / (L_m * k)) / math.log(term_inside_ln)
        T_corrected_rounded = int(round(T_corrected))

        print("The provided measurements are inconsistent with Planck's Law for black-body radiation.")
        print("The measured temperature of 9000 K is unusually low for a DB-class white dwarf.")
        print("Assuming the spectral radiance and wavelength measurements are correct, the expected temperature can be calculated.")
        print("\nThe equation to find the corrected temperature (T_corr) is:")
        print("T_corr = (h * c) / (λ * k * ln((2 * h * c^2) / (λ^5 * B) + 1))")
        print("\nPlugging in the values:")
        print(f"T_corr = ({h:.4e} * {c:.4e}) / ({L_m:.4e} * {k:.4e} * ln((2 * {h:.4e} * ({c:.4e})^2) / (({L_m:.4e})^5 * {B_measured:.4e}) + 1))")
        print(f"\nCalculated Corrected Temperature: {T_corrected:.2f} K")
        print(f"Rounded Corrected Temperature: {T_corrected_rounded} K")
        
        # Final answer in the format Xv
        print(f"\nT{T_corrected_rounded}")

solve_pandora_problem()