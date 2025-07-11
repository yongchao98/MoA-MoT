import math

def calculate_dynamic_structure_factor_denominator():
    """
    Calculates the denominator of the dynamic structure factor S(k, ω)
    for a given set of nuclear parameters and transform variables.
    """
    # --- System Parameters for a hypothetical subcritical system ---

    # Reactivity (ρ), a measure of criticality. Negative for subcritical.
    rho = -0.01

    # Effective delayed neutron fraction (β)
    beta = 0.0065

    # Prompt neutron generation time (Λ) in seconds
    Lambda = 2.0e-5

    # Neutron diffusion coefficient (D) in cm²/s
    D = 2.0e5

    # --- Fourier Transform Variables ---

    # Wave number (k) in cm⁻¹
    k = 0.1

    # Frequency (f) in Hz
    f = 100.0

    # --- Calculations ---

    # Calculate Rossi-alpha (α), the prompt neutron decay constant in s⁻¹
    alpha = (beta - rho) / Lambda

    # Convert frequency in Hz to angular frequency (ω) in rad/s
    omega = 2 * math.pi * f

    # Calculate the components of the denominator
    alpha_term = alpha
    Dk2_term = D * k**2
    omega_term = omega
    
    numerator_of_denominator = alpha_term + Dk2_term
    denominator = numerator_of_denominator**2 + omega_term**2

    # --- Output ---

    print("The dynamic structure factor S(k, ω) is given by S(k, ω) ∝ 1 / [ (α + D*k²)² + ω² ]\n")
    
    print("System Parameters:")
    print(f"  Reactivity (ρ): {rho}")
    print(f"  Delayed Neutron Fraction (β): {beta}")
    print(f"  Generation Time (Λ): {Lambda:.1e} s")
    print(f"  Diffusion Coefficient (D): {D:.1e} cm²/s")
    print(f"  => Calculated Rossi-alpha (α): {alpha:.2f} s⁻¹\n")

    print("Fourier Variables:")
    print(f"  Wave number (k): {k} cm⁻¹")
    print(f"  Frequency (f): {f} Hz")
    print(f"  => Angular Frequency (ω): {omega:.2f} rad/s\n")
    
    print("Final Equation Calculation:")
    print("S(k, ω) ∝ 1 / [ (α + D*k²)² + ω² ]")
    # Outputting each number in the final equation
    print(f"S(k={k}, ω={omega:.2f}) ∝ 1 / [ ({alpha_term:.2f} + {D:.1e}*{k}² )² + {omega_term:.2f}² ]")
    print(f"S(k={k}, ω={omega:.2f}) ∝ 1 / [ ({alpha_term:.2f} + {Dk2_term:.2f})² + {omega_term:.2f}² ]")
    print(f"S(k={k}, ω={omega:.2f}) ∝ 1 / [ ({numerator_of_denominator:.2f})² + {omega_term:.2f}² ]")
    print(f"S(k={k}, ω={omega:.2f}) ∝ 1 / [ {numerator_of_denominator**2:.2e} + {omega_term**2:.2e} ]")
    print(f"S(k={k}, ω={omega:.2f}) ∝ 1 / {denominator:.2e}")

if __name__ == "__main__":
    calculate_dynamic_structure_factor_denominator()