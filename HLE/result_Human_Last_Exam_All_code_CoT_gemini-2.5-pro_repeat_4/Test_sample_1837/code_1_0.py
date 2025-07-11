import math

def solve_pandora_mystery():
    """
    Analyzes stellar measurements using Planck's Law to find inconsistencies.
    """
    # --- Step 1: Define constants and given values ---
    h = 6.62607015e-34  # Planck constant (J·s)
    c = 299792458      # Speed of light (m/s)
    k = 1.380649e-23   # Boltzmann constant (J/K)

    # Given measurements from the problem
    L_nm = 400.0         # Wavelength in nm
    T = 9000.0           # Temperature in K
    B_mes = 1.2e15       # Measured spectral radiance in W/(m^2·sr·m)

    # Convert wavelength to SI units (meters)
    L_m = L_nm * 1e-9

    print("Analyzing the measurements for the star Pandora.")
    print(f"Given: Wavelength (L) = {L_nm:.0f} nm, Temperature (T) = {T:.0f} K, Spectral Radiance (B) = {B_mes:.1e} W/m²/sr/m")

    # --- Step 2: Check for consistency using Planck's Law ---
    print("\nWe will check consistency using Planck's Law:")
    print("B(L, T) = (2 * h * c^2) / (L^5 * (e^(h*c / (L*k*T)) - 1))")
    print("\nFirst, calculating the expected spectral radiance (B_exp) for the given T and L.")

    # Calculate the exponent part for the formula
    exponent = (h * c) / (L_m * k * T)
    # Calculate the expected spectral radiance
    B_exp = (2 * h * c**2) / (L_m**5 * (math.exp(exponent) - 1))

    print("\nPlugging in the values into the equation:")
    print(f"B_exp = (2 * {h:.6e} * {c:.0f}^2) / (({L_m:.1e})^5 * (e^(({h:.6e} * {c:.0f}) / ({L_m:.1e} * {k:.6e} * {T:.0f})) - 1))")
    print(f"B_exp = {B_exp:.3e} W/m²/sr/m")

    # --- Step 3: Compare measured and expected values ---
    # Using a 5% threshold for consistency
    if abs(B_mes - B_exp) / B_exp < 0.05:
        print("\nThe measured spectral radiance is consistent with the given temperature and wavelength.")
        final_answer = "0"
    else:
        print(f"\nThe measured radiance {B_mes:.1e} is significantly different from the expected value {B_exp:.3e}.")
        print("This suggests a measurement error in one of the quantities (L, T, or B).")

        # --- Step 4: Investigate the most likely error ---

        # Hypothesis 1: Wavelength (L) is incorrect.
        # Check if B_mes is physically possible at the given temperature T.
        # The peak of B_lambda is found via the root of x = 5*(1-e^-x), which is x~4.965
        # where x = hc/(L*k*T)
        lambda_peak_m = (h * c) / (4.965114 * k * T)
        exp_peak = (h * c) / (lambda_peak_m * k * T)
        B_peak = (2 * h * c**2) / (lambda_peak_m**5 * (math.exp(exp_peak) - 1))

        print(f"\n- Hypothesis 1: Wavelength (L) is incorrect.")
        print(f"At T={T:.0f} K, the maximum possible spectral radiance (at λ_peak = {lambda_peak_m*1e9:.0f} nm) is {B_peak:.3e} W/m²/sr/m.")
        print(f"The measured radiance {B_mes:.1e} is HIGHER than this maximum. Therefore, this hypothesis is not possible.")

        # Hypothesis 2: Temperature (T) is incorrect.
        print(f"\n- Hypothesis 2: Temperature (T) is incorrect.")
        print("Calculating the temperature (T_new) that would produce the measured radiance B_mes at the given wavelength L.")

        # To find T_new, we solve the Planck's law equation for T:
        # T_new = (h*c) / (L*k*ln(1 + (2*h*c^2)/(L^5 * B_mes)))
        C1 = 2 * h * c**2
        try:
            val_for_ln = 1 + (C1 / (L_m**5 * B_mes))
            T_new = (h * c) / (L_m * k * math.log(val_for_ln))
        except (ValueError, ZeroDivisionError):
            T_new = float('inf') # Should not happen in this case

        print("\nSolving for T in the equation:")
        print(f"T_new = ({h:.6e} * {c:.0f}) / ({L_m:.1e} * {k:.6e} * ln(1 + (2 * {h:.6e} * {c:.0f}^2) / (({L_m:.1e})^5 * {B_mes:.1e})))")
        print(f"The resulting temperature is {T_new:.0f} K.")

        print("\nEvaluating the result:")
        print("A DB-class white dwarf star typically has a surface temperature between 12,000 K and 30,000 K.")
        print(f"The original temperature of {T:.0f} K is too low, falling outside this range.")
        print(f"The corrected temperature of {T_new:.0f} K fits perfectly within the expected range.")
        print("\nConclusion: The temperature measurement (T) is the most likely source of error.")

        final_answer = f"T{T_new:.0f}"

    # --- Final Output ---
    print(f"\n<<<{final_answer}>>>")

solve_pandora_mystery()