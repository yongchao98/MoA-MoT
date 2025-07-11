import math

def solve_exponent_order():
    """
    Calculates the order of the first non-vanishing contribution to the
    critical exponent nu in the perturbative expansion of phi^4 theory.
    """

    print("This program determines the order of the coupling constant 'u' for the first correction to the critical exponent nu.")
    print("-" * 80)
    
    print("\nStep 1: The RG relation for the critical exponent nu")
    print("In phi^4 theory, the exponent nu(u) is given by:")
    print("  1 / nu(u) = 2 - gamma(u)")
    print("where gamma(u) is the anomalous dimension of the phi^2 operator.")

    print("\nStep 2: The one-loop perturbative expansion for gamma(u)")
    print("The one-loop calculation for gamma(u) gives it as a power series in u, starting with a linear term:")
    # For a single-component scalar field (N=1), C1 = 1 / (8 * pi^2).
    # We will use this value for the demonstration.
    C1 = 1 / (8 * math.pi**2)
    gamma_coeffs = [0, C1]  # Represents gamma(u) = 0*u^0 + C1*u^1 + ...
    print(f"  gamma(u) = {gamma_coeffs[0]:.1f} + {gamma_coeffs[1]:.5f}*u + O(u^2)")

    print("\nStep 3: Calculating the power series for 1 / nu(u)")
    # From the relation in Step 1: 1/nu(u) = 2 - (C1*u + ...)
    inv_nu_coeffs = [2.0 - gamma_coeffs[0], -gamma_coeffs[1]]
    print(f"  1 / nu(u) = {inv_nu_coeffs[0]:.1f} + ({inv_nu_coeffs[1]:.5f})*u + O(u^2)")

    print("\nStep 4: Calculating the power series for nu(u) via series inversion")
    print("To find nu(u), we invert the series for 1 / nu(u).")
    # If a function is f(u) = a0 + a1*u, its inverse is 1/f(u) = (1/a0) - (a1/a0^2)*u + ...
    a0 = inv_nu_coeffs[0]
    a1 = inv_nu_coeffs[1]
    nu_c0 = 1 / a0
    nu_c1 = -a1 / (a0**2)
    nu_coeffs = [nu_c0, nu_c1]
    print(f"The expansion for nu(u) is: nu(u) = {nu_coeffs[0]:.3f} + {nu_coeffs[1]:.5f}*u + O(u^2)")
    
    print("\nStep 5: Identifying the order of the first correction")
    print(f"The mean-field value (at u=0) is nu(0) = {nu_coeffs[0]:.3f}")
    
    final_equation_nu0 = nu_coeffs[0]
    final_equation_nu1 = nu_coeffs[1]
    order = 1
    
    print("The first correction to this value is the next term in the series:")
    print(f"Correction term = ({final_equation_nu1:.5f}) * u^{order}")
    print("\nThe order of the coupling 'u' in this initial non-vanishing contribution is the exponent of u.")
    
    print("\n--- Final Equation and Result ---")
    print("The first-order expansion for nu(u) is: ")
    print(f"nu(u) ≈ {final_equation_nu0:.3f} + ({final_equation_nu1:.5f}) * u^{order}")
    print("\nThe requested order is the exponent of u in the correction term.")
    print(f"Order = {order}")


solve_exponent_order()