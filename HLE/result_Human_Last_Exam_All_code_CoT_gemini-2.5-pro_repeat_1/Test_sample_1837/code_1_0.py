import math

def solve_star_measurement():
    """
    Validates astronomical measurements for a star using Planck's Law.
    
    Given:
    - Wavelength (L): 400 nm
    - Spectral Radiance (B): 1.2e15 W·m⁻²·sr⁻¹·m⁻¹
    - Surface Temperature (T): 9000 K
    
    This function checks for inconsistencies and identifies the most likely
    erroneous measurement and its corrected value.
    """
    
    # --- Constants ---
    h = 6.62607015e-34  # Planck constant (J·s)
    c = 299792458      # Speed of light (m/s)
    k = 1.380649e-23   # Boltzmann constant (J/K)

    # --- Given Measurements ---
    lambda_nm = 400.0
    lambda_m = lambda_nm * 1e-9  # Convert wavelength to meters
    B_measured = 1.2e15
    T_given = 9000.0

    # --- Analysis ---
    # First, calculate the theoretical maximum spectral radiance for the given temperature.
    # This occurs at the peak of the black-body curve, given by Wien's Law.
    b_wien = 2.898e-3  # Wien's displacement constant in m·K
    lambda_peak_m = b_wien / T_given
    
    # Calculate max spectral radiance using Planck's law at this peak wavelength
    exp_term_max = math.exp((h * c) / (lambda_peak_m * k * T_given))
    B_max_possible = (2 * h * c**2) / (lambda_peak_m**5 * (exp_term_max - 1))

    # The measured radiance (1.2e15) is significantly higher than the maximum possible
    # radiance at 9000 K (~2.4e14). This indicates the temperature measurement is incorrect.
    # We will now calculate the correct temperature assuming B and L are correct.

    # --- Recalculate Temperature ---
    # Rearrange Planck's Law to solve for T:
    # T = (h*c) / (λ*k*ln( (2*h*c^2)/(B*λ^5) + 1 ))
    
    # Break down the calculation for clarity
    numerator_main = h * c
    denominator_in_ln = B_measured * (lambda_m**5)
    numerator_in_ln = 2 * h * c**2
    term_in_ln = numerator_in_ln / denominator_in_ln
    ln_result = math.log(term_in_ln + 1)
    
    T_corrected = numerator_main / (lambda_m * k * ln_result)
    T_corrected_rounded = round(T_corrected)

    # --- Output Results ---
    print("Recalculating Temperature (T) assuming Spectral Radiance (B) and Wavelength (λ) are correct.")
    print("Using Planck's Law: T = (h * c) / (λ * k * ln( (2 * h * c^2 / (B * λ^5)) + 1 ))")
    print("\n--- Calculation Breakdown ---")

    print(f"h (Planck constant) = {h:.4e} J·s")
    print(f"c (Speed of light) = {c:.4e} m/s")
    print(f"k (Boltzmann constant) = {k:.4e} J/K")
    print(f"λ (Wavelength) = {lambda_m:.4e} m")
    print(f"B (Spectral Radiance) = {B_measured:.1e} W·m⁻²·sr⁻¹·m⁻¹")

    print("\n--- Equation with values ---")
    print(f"Numerator (h * c) = {h * c:.4e}")
    print(f"Term in ln (2*h*c^2 / (B*λ^5)) = ({numerator_in_ln:.4e}) / ({denominator_in_ln:.4e}) = {term_in_ln:.4f}")
    print(f"ln(...) = ln({term_in_ln:.4f} + 1) = {ln_result:.4f}")
    
    print(f"\nT = ({h * c:.4e}) / ({lambda_m:.4e} * {k:.4e} * {ln_result:.4f})")
    print(f"T = {T_corrected:.0f} K")

    # --- Final Answer ---
    final_answer = f"T{T_corrected_rounded}"
    print(f"\n<<<{final_answer}>>>")

# Run the solution
solve_star_measurement()