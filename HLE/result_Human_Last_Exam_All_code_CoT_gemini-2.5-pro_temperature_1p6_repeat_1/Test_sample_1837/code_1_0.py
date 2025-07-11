import math

def solve_pandora_problem():
    """
    Analyzes the consistency of stellar measurements using Planck's Law and
    calculates the corrected temperature.
    """
    # --- Physical constants ---
    h = 6.62607015e-34  # Planck's constant (J*s)
    c = 299792458      # Speed of light (m/s)
    k_B = 1.380649e-23   # Boltzmann constant (J/K)

    # --- Given measured values ---
    lambda_val = 400e-9      # Wavelength in meters (400 nm)
    T_measured = 9000.0      # Surface temperature in Kelvin
    B_measured = 1.2e15      # Spectral radiance in W/(m^2*sr*m)

    # --- Check for consistency by calculating what B should be for the given T and λ ---
    # B(λ, T) = (2 * h * c^2 / λ^5) * (1 / (e^(h*c / (λ*k_B*T)) - 1))
    try:
        exponent_val = (h * c) / (lambda_val * k_B * T_measured)
        B_theoretical = (2 * h * c**2) / (lambda_val**5 * (math.exp(exponent_val) - 1))
    except OverflowError:
        B_theoretical = 0.0
    
    print(f"Checking consistency: For T={T_measured} K and λ={lambda_val*1e9:.0f} nm, the theoretical spectral radiance is {B_theoretical:.2e} W/(m^2*sr*m).")
    print(f"This is significantly different from the measured value of {B_measured:.2e} W/(m^2*sr*m).\n")
    
    print("The measured temperature of 9000 K is outside the typical range for a DB-class white dwarf (12,000 K - 30,000 K).")
    print("Therefore, temperature (T) is the quantity most likely in error.\n")

    # --- Calculate the corrected temperature (T) assuming B and λ are correct ---
    print("Calculating the expected correct temperature assuming B and λ are correct.")
    print("Using the formula: T = (h*c) / [λ*k_B * ln(1 + (2*h*c^2)/(λ^5*B))]\n")
    
    # Pre-calculating parts of the equation for clarity in the printout
    hc = h * c
    two_h_c_sq = 2 * h * c**2
    lambda_5_B = lambda_val**5 * B_measured
    lambda_k = lambda_val * k_B
    term_in_ln = 1 + (two_h_c_sq / lambda_5_B)
    
    # Perform the final calculation
    T_corrected = hc / (lambda_k * math.log(term_in_ln))
    T_corrected_rounded = round(T_corrected)

    # Print the equation with the values plugged in
    print("Substituting the measured values into the equation:")
    print(f"h*c = {h:.4e} * {c:.4e} = {hc:.4e}")
    print(f"2*h*c^2 = 2 * {h:.4e} * ({c:.4e})^2 = {two_h_c_sq:.4e}")
    print(f"λ^5 * B = ({lambda_val:.4e})^5 * {B_measured:.4e} = {lambda_5_B:.4e}")
    print(f"λ*k_B = {lambda_val:.4e} * {k_B:.4e} = {lambda_k:.4e}")
    print(f"\nT = {hc:.4e} / [{lambda_k:.4e} * ln(1 + {two_h_c_sq:.4e} / {lambda_5_B:.4e})]")
    print(f"T = {hc:.4e} / [{lambda_k:.4e} * ln({term_in_ln:.4f})]")
    print(f"T = {hc:.4e} / [{lambda_k:.4e} * {math.log(term_in_ln):.4f}]")
    print(f"T = {T_corrected:.2f} K\n")

    print(f"The expected correct value for T, rounded to the nearest unit, is: {T_corrected_rounded} K")

solve_pandora_problem()