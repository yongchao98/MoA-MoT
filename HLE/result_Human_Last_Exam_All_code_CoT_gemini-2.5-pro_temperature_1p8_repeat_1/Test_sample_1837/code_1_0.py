import math

def solve_pandora_problem():
    """
    This function verifies the measurements of the star Pandora and finds the corrected value
    for the most likely erroneous measurement.
    """

    # Physical constants in SI units
    h = 6.62607015e-34  # Planck constant (J·s)
    c = 299792458      # Speed of light (m/s)
    k = 1.380649e-23   # Boltzmann constant (J/K)

    # Given measurements
    lam = 400e-9         # Wavelength in meters (from 400 nm)
    B = 1.2e15           # Measured spectral radiance (W·m⁻²·sr⁻¹·m⁻¹)
    
    # As reasoned in the plan, temperature (T) is the most likely source of error.
    # We solve for T by inverting Planck's Law:
    # T = (h * c) / (λ * k * ln((2 * h * c^2) / (λ^5 * B) + 1))
    
    # Let's output the numbers involved in the calculation as requested.
    # The equation for the corrected temperature T is derived from B(λ,T).
    # We substitute the following numbers:
    
    hc = h * c
    two_hc_sq = 2 * h * c**2
    lam_5_B = (lam**5) * B
    
    term_in_log = two_hc_sq / lam_5_B + 1
    
    lam_k_ln = lam * k * math.log(term_in_log)

    T_corrected = hc / lam_k_ln
    T_final = round(T_corrected)

    print("Equation for corrected temperature: T = (h*c) / (λ*k*ln( (2*h*c²)/(λ⁵*B) + 1 ))")
    print(f"h*c = {hc:.4g}")
    print(f"2*h*c² = {two_hc_sq:.4g}")
    print(f"λ = {lam:.4g}")
    print(f"k = {k:.4g}")
    print(f"B = {B:.4g}")
    print(f"λ⁵*B = {lam_5_B:.4g}")
    print(f"T = {hc:.4g} / ( {lam:.4g} * {k:.4g} * ln({term_in_log:.4f}) ) = {T_final}")
    
    # Final answer in the required format Xv
    print(f"T{T_final}")

solve_pandora_problem()