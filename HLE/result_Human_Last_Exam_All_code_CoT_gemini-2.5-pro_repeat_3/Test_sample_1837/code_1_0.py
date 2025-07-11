import math

def solve_star_measurement():
    """
    Analyzes the measurements of the white dwarf star Pandora and identifies the likely error.
    """
    # Constants
    h = 6.62607015e-34  # Planck's constant (J·s)
    c = 2.99792458e8    # Speed of light (m/s)
    k_B = 1.380649e-23   # Boltzmann constant (J/K)

    # Given measurements
    lam_nm = 400
    lam_m = lam_nm * 1e-9  # Wavelength in meters
    T_measured = 9000      # Temperature in Kelvin
    B_measured = 1.2e15    # Spectral radiance in W/(m^2·sr·m)

    print("Step 1: Verify the given measurements using Planck's Law.")
    print("Planck's Law: B(λ, T) = (2 * h * c^2) / (λ^5 * (exp((h * c) / (λ * k_B * T)) - 1))")
    print("-" * 50)
    
    # --- Calculate expected Spectral Radiance (B) for the given T and L ---
    print("Calculating the expected spectral radiance (B_expected) for T = 9000 K and L = 400 nm.")
    
    exponent_val = (h * c) / (lam_m * k_B * T_measured)
    denominator = (lam_m**5) * (math.exp(exponent_val) - 1)
    B_expected = (2 * h * c**2) / denominator

    # Print the equation with numbers
    print("\nEquation for B_expected:")
    print(f"B = (2 * {h:.3e} * ({c:.3e})^2) / (({lam_m:.3e})^5 * (exp(({h:.3e} * {c:.3e}) / ({lam_m:.3e} * {k_B:.3e} * {T_measured})) - 1))")
    print(f"B_expected = {B_expected:.3e} W/m²/sr/m")
    print(f"B_measured = {B_measured:.3e} W/m²/sr/m")

    # --- Compare and analyze ---
    print("\nStep 2: Compare B_expected with B_measured and analyze the discrepancy.")
    print("-" * 50)
    print(f"The calculated radiance ({B_expected:.3e}) and the measured radiance ({B_measured:.3e}) differ significantly.")
    print("This indicates an error in at least one of the initial measurements (L, T, or B).")
    
    # --- Check for plausibility ---
    # The peak radiance for a 9000K black body is lower than the measured radiance,
    # which means L cannot be the only error. The error must be in T or B.
    # The crucial clue is the star's classification.
    print("\nAstronomical context: DB-class white dwarfs typically have temperatures from 12,000 K to 30,000 K.")
    print(f"The given temperature of {T_measured} K is unusually low for this star type.")

    # --- Assume L and B are correct, and calculate the correct Temperature (T) ---
    print("\nStep 3: Calculate the required temperature (T_corrected) to match the measured radiance.")
    print("-" * 50)
    
    # Rearranging Planck's law to solve for T:
    # T = (h * c) / (λ * k_B * ln(1 + (2 * h * c^2) / (λ^5 * B)))
    try:
        log_term_arg = 1 + (2 * h * c**2) / (lam_m**5 * B_measured)
        if log_term_arg <= 0:
            raise ValueError("Invalid argument for logarithm, radiance value is too high for any real temperature.")
        
        T_corrected = (h * c) / (lam_m * k_B * math.log(log_term_arg))

        print("Equation for T_corrected:")
        print(f"T = ({h:.3e} * {c:.3e}) / ({lam_m:.3e} * {k_B:.3e} * ln(1 + (2 * {h:.3e} * ({c:.3e})^2) / (({lam_m:.3e})^5 * {B_measured:.3e})))")
        print(f"T_corrected = {T_corrected:.0f} K")
        
        print(f"\nThis calculated temperature of {T_corrected:.0f} K falls well within the typical range for a DB-class white dwarf.")
        print("Therefore, the temperature (T) is the most likely quantity to have a measurement error.")
        
        final_answer = f"T{T_corrected:.0f}"
        print("\nFinal Answer:")
        print(final_answer)

    except ValueError as e:
        print(f"Error during calculation: {e}")

solve_star_measurement()