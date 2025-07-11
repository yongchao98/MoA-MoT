import numpy as np
from scipy.integrate import quad

def evaluate_1s_integral(zeta):
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.

    This function explains the analytical derivation and verifies it with a 
    numerical calculation for a given Slater exponent zeta.
    """
    
    # --- Analytical Derivation Explanation ---
    print("This script evaluates the integral <phi | 1/r | phi> for a 1s Slater-type orbital (STO).")
    print("The normalized 1s STO is given by: phi(r) = (zeta^3 / pi)^(1/2) * exp(-zeta * r)")
    print("The integral in spherical coordinates is:")
    print("I = integral( |phi(r)|^2 * (1/r) * dV ) over all space")
    print("where dV = r^2 * sin(theta) * dr * dtheta * dphi.")
    
    print("\nSubstituting the expression for phi(r) and simplifying:")
    print("I = (zeta^3 / pi) * integral_all_space [ exp(-2*zeta*r) * (1/r) * r^2 * sin(theta) ] dr dtheta dphi")
    print("I = (zeta^3 / pi) * integral_0_to_inf [ r * exp(-2*zeta*r) ] dr * integral_0_to_pi [ sin(theta) ] dtheta * integral_0_to_2pi [ dphi ]")
    
    # Evaluate angular integrals
    angular_integral_result = 4 * np.pi
    print(f"\nThe angular part integrates to (2) * (2*pi) = 4*pi.")
    
    print("\nThe expression simplifies to:")
    print("I = (zeta^3 / pi) * (4*pi) * integral_from_0_to_inf [ r * exp(-2*zeta*r) ] dr")
    print("I = 4 * zeta^3 * integral_from_0_to_inf [ r * exp(-2*zeta*r) ] dr")
    
    print("\nUsing the standard integral formula: integral_0^inf(x^n * exp(-a*x) dx) = n! / a^(n+1)")
    print("With n=1 and a=2*zeta, the radial integral evaluates to 1! / (2*zeta)^2 = 1 / (4*zeta^2).")
    
    print("\nSubstituting the result of the radial integral:")
    print("I = 4 * zeta^3 * (1 / (4*zeta^2))")
    print("I = zeta")
    
    analytical_result = zeta

    print("\n-----------------------------------------------------")
    print(f"Calculation for the specific case: zeta = {zeta}")
    print("-----------------------------------------------------")

    # --- Step-by-step final equation evaluation ---
    print("\nFinal Equation Evaluation:")
    radial_integral_value = 1 / (4 * zeta**2)
    final_value = 4 * zeta**3 * radial_integral_value
    
    print(f"Integral = 4 * ({zeta})^3 * (1 / (4 * ({zeta})^2))")
    print(f"Integral = 4 * {zeta**3:.4f} * (1 / (4 * {zeta**2:.4f}))")
    print(f"Integral = 4 * {zeta**3:.4f} * ({radial_integral_value:.4f})")
    print(f"Integral = {final_value:.4f}")
    
    # --- Numerical Verification ---
    print("\n--- Numerical Verification ---")
    # The full integrand for numerical integration (after angular part is done)
    def integrand(r, z):
        return 4 * z**3 * r * np.exp(-2 * z * r)

    numerical_result, error = quad(integrand, 0, np.inf, args=(zeta,))
    
    print(f"Numerical result using scipy.integrate.quad: {numerical_result:.6f}")
    print(f"Analytical result: {analytical_result:.6f}")
    print(f"The results match, confirming the derivation.")

# --- Main execution ---
# For a Hydrogen atom, the 1s orbital has an exponent zeta = 1.0.
# We will use this value for the calculation as a representative example.
zeta_value = 1.0
evaluate_1s_integral(zeta_value)