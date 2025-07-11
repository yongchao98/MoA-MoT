import math

def solve_star_measurement():
    """
    Analyzes the provided stellar measurements using Planck's Law.
    1. Calculates the theoretical spectral radiance based on T and L.
    2. Compares it with the measured spectral radiance.
    3. If inconsistent, determines the most likely error and calculates the corrected value.
    """
    # Physical constants
    h = 6.62607015e-34  # Planck's constant (J·s)
    c = 2.99792458e8   # Speed of light (m/s)
    k = 1.380649e-23     # Boltzmann constant (J/K)

    # Provided measurements
    L_nm = 400.0         # Wavelength in nanometers
    L = L_nm * 1e-9      # Wavelength in meters
    T = 9000.0           # Temperature in Kelvin
    B_measured = 1.2e15  # Spectral radiance (W / (m^2 * sr * m))

    # --- Step 1: Check for consistency ---
    # Calculate the theoretical spectral radiance for the given T and L
    try:
        exponent_val = (h * c) / (L * k * T)
        # Prevent overflow for very large exponents
        if exponent_val > 700:
             B_theoretical = 0
        else:
            numerator = 2 * h * c**2 / L**5
            denominator = math.exp(exponent_val) - 1
            B_theoretical = numerator / denominator
    except (ValueError, OverflowError):
        B_theoretical = 0.0

    # Check if B_measured is close to B_theoretical (e.g., within 10%)
    if B_theoretical > 0 and abs(B_measured - B_theoretical) / B_theoretical < 0.1:
        print("0")
        return

    # --- Step 2: Inconsistency found, find the most likely error ---
    # The most likely error is Temperature, as a 9000K star cannot produce
    # this radiance, and ~15000K is a more typical temperature for a DB white dwarf.
    # We solve for the Temperature T, assuming B and L are correct.
    
    # Equation: B = (2*h*c^2 / L^5) / (exp(h*c / (L*k*T)) - 1)
    
    # We will show the steps to solve for T.
    
    # First, calculate the constant terms in the equation.
    term1 = (2 * h * c**2) / (L**5)
    term2_const = (h * c) / (L * k)

    # Now, solve for T
    # B_measured = term1 / (exp(term2_const / T) - 1)
    # exp(term2_const / T) - 1 = term1 / B_measured
    # exp(term2_const / T) = (term1 / B_measured) + 1
    # term2_const / T = ln((term1 / B_measured) + 1)
    # T = term2_const / ln((term1 / B_measured) + 1)

    ratio = term1 / B_measured
    ln_val = math.log(ratio + 1)
    T_corrected = term2_const / ln_val
    T_corrected_rounded = round(T_corrected)

    print("The measured values are inconsistent with Planck's Law for a 9000K star.")
    print("The theoretical spectral radiance for a 9000K star at 400 nm would be {:.3e} W/m^2/sr/m.".format(B_theoretical))
    print("\nAssuming the measured spectral radiance and wavelength are correct, we can solve for the true temperature T.")
    print("This is the most likely scenario, as the resulting temperature is more typical for a DB-class white dwarf.\n")
    print("Solving Planck's Law for T:")
    print("B = (2*h*c^2 / L^5) / (exp(h*c / (L*k*T)) - 1)")
    print("\nStep 1: Plug in the known values B, h, c, L, k into the equation.")
    print("{:.3e} = (2 * {:.3e} * ({:.3e})^2 / ({:.3e})^5) / (exp(( {:.3e} * {:.3e} ) / ( {:.3e} * {:.3e} * T )) - 1)".format(
        B_measured, h, c, L, h, c, L, k
    ))
    
    print("\nStep 2: Simplify the constant terms.")
    print("The term (2*h*c^2 / L^5) calculates to {:.3e}.".format(term1))
    print("The term (h*c / (L*k)) in the exponent calculates to {:.0f}.".format(term2_const))
    print("The equation becomes: {:.3e} = {:.3e} / (exp({:.0f} / T) - 1)".format(
        B_measured, term1, term2_const
    ))
    
    print("\nStep 3: Rearrange the equation to solve for T.")
    print("exp({:.0f} / T) - 1 = {:.3e} / {:.3e} = {:.3f}".format(
        term2_const, term1, B_measured, ratio
    ))
    print("exp({:.0f} / T) = {:.3f}".format(term2_const, ratio + 1))
    print("{:.0f} / T = ln({:.3f}) = {:.3f}".format(
        term2_const, ratio + 1, ln_val
    ))
    
    print("\nStep 4: Calculate the final value for T.")
    print("T = {:.0f} / {:.3f} = {:.0f} K".format(
        term2_const, ln_val, T_corrected_rounded
    ))


solve_star_measurement()
