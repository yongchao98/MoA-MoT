import math

def solve_pandora_problem():
    """
    This script checks the consistency of stellar measurements using Planck's Law
    and identifies the most likely erroneous value.
    """
    # Pandora's data from the problem description
    # Wavelength in nm
    L_nm = 400.0
    # Spectral radiance in W/(m^2*sr*m)
    B_meas = 1.2e15
    # Estimated Temperature in K
    T_est = 9000.0

    # Physical constants
    h = 6.62607015e-34  # Planck's constant (J*s)
    c = 299792458      # Speed of light (m/s)
    k = 1.380649e-23   # Boltzmann constant (J/K)

    # Convert wavelength from nanometers to meters for calculations
    L_m = L_nm * 1e-9

    # Based on initial analysis, the measured spectral radiance (B_meas) is
    # physically impossible at the estimated temperature (T_est).
    # Therefore, the temperature is the most likely quantity with an error.
    # We will calculate the correct temperature (T') assuming B_meas and L_m are accurate.
    # The formula, derived from Planck's Law, is:
    # T' = (h * c) / (λ * k * ln( (2 * h * c^2) / (B * λ^5) + 1 ))

    # Calculate the individual terms for clarity
    term_in_ln_numerator = 2 * h * c**2
    term_in_ln_denominator = B_meas * (L_m**5)

    # A potential issue is if the argument to log is negative, but B is always positive
    # and the term (2*h*c^2)/(B*L^5) will be positive.
    term_in_ln = term_in_ln_numerator / term_in_ln_denominator
    log_term = math.log(term_in_ln + 1)
    
    numerator = h * c
    denominator = L_m * k * log_term

    # Calculate the corrected temperature
    T_corr = numerator / denominator
    T_corr_rounded = round(T_corr)

    # As requested, output the numbers in the final equation
    print("The corrected temperature T' is found using the following equation derived from Planck's Law:")
    print("T' = (h * c) / (λ * k * ln( (2 * h * c^2) / (B * λ^5) + 1 ))")
    
    print("\nPlugging in the measured values and physical constants:")
    print(f"h (Planck's constant) = {h} J*s")
    print(f"c (Speed of light) = {c} m/s")
    print(f"k (Boltzmann constant) = {k} J/K")
    print(f"λ (Wavelength) = {L_m} m")
    print(f"B (Spectral Radiance) = {B_meas} W/(m^2*sr*m)")

    print("\nThe equation with the numbers is:")
    print(f"T' = ({numerator:.6e}) / ({L_m} * {k:.6e} * ln( ({term_in_ln_numerator:.6e}) / ({B_meas} * ({L_m})**5) + 1 ))")
    print(f"T' = {T_corr:.2f} K")

    print(f"\nSince the calculated temperature ({T_corr_rounded} K) is significantly different from the estimate ({T_est} K), the temperature is the most likely erroneous quantity.")
    
    # The final answer must be in the format Xv
    print(f"\nFinal Answer: T{T_corr_rounded}")

solve_pandora_problem()