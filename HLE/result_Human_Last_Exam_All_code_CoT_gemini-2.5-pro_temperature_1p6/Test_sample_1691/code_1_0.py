import numpy as np
from scipy.special import beta

def solve_integral_approximation():
    """
    Calculates the coefficients for the asymptotic expansion of the integral I(epsilon).
    """
    # Problem parameters from the denominator:
    # f(x) = c*x^n + d*x^m + e*x^p
    # c=9.0, n=5.0
    # d=5.0, m=6.0
    # e=9.0, p=8.0
    # Integration limit
    # L=15.0

    c, n = 9.0, 5.0
    d, m = 5.0, 6.0
    e, p_power = 9.0, 8.0
    L = 15.0

    # --- Coefficient C1 for the epsilon^(-4/5) term ---
    # The integral is of the form integral(dy / (1 + c*y^n))
    # After substitution u^n = c*y^n, this becomes c^(-1/n) * integral(du / (1+u^n))
    # integral(du / (1+u^n)) from 0 to inf is pi / (n * sin(pi/n))
    integral_part_C1 = np.pi / (n * np.sin(np.pi / n))
    C1 = (c**(-1/n)) * integral_part_C1

    # --- Coefficient C2 for the epsilon^(-3/5) term ---
    # This term comes from the -d*x^m term in the expansion
    # The integral is -d * integral(y^m / (1+c*y^n)^2) dy
    # This can be evaluated using the Beta function B(x, y)
    # The integral is (1/n)*c^-(m+1)/n * B((m+1)/n, 2-(m+1)/n)
    x_beta2 = (m + 1) / n
    y_beta2 = 2 - x_beta2
    integral_part_C2 = (1 / n) * (c**(-(m + 1) / n)) * beta(x_beta2, y_beta2)
    C2 = -d * integral_part_C2

    # --- Coefficient C3 for the epsilon^(-1/5) term ---
    # This term comes from the -e*x^p term in the expansion
    x_beta3 = (p_power + 1) / n
    y_beta3 = 2 - x_beta3
    integral_part_C3 = (1 / n) * (c**(-(p_power + 1) / n)) * beta(x_beta3, y_beta3)
    C3 = -e * integral_part_C3

    # --- Constant term C0 ---
    # This arises from correcting the integrals for the finite limit L
    C0_term1 = -1 / (c * (n - 1) * (L**(n - 1)))
    C0_term2 = d / (c**2 * (2 * n - m - 1) * L**(2 * n - m - 1))
    C0_term3 = e / (c**2 * (2 * n - p_power - 1) * L**(2 * n - p_power - 1))
    C0 = C0_term1 + C0_term2 + C0_term3

    print("The analytical formula for the integral I(epsilon) for small epsilon is an asymptotic expansion:")
    print("I(epsilon) = C1 * epsilon^(-4/5) + C2 * epsilon^(-3/5) + C3 * epsilon^(-1/5) + C0 + O(epsilon^(1/5))")
    print("\nCalculated coefficients:")
    print(f"C1 (for epsilon^-0.8) = {C1:.6f}")
    print(f"C2 (for epsilon^-0.6) = {C2:.6f}")
    print(f"C3 (for epsilon^-0.2) = {C3:.6f}")
    print(f"C0 (constant term)  = {C0:.6f}")
    
    print("\nFinal formula:")
    print(f"I(epsilon) approx = {C1:.4f} * epsilon^(-4/5) + ({C2:.4f}) * epsilon^(-3/5) + ({C3:.4f}) * epsilon^(-1/5) + {C0:.4f}")
    
    print("\nBreakdown of the constant term C0:")
    print(f"Contribution from the leading term = {C0_term1:.6f}")
    print(f"Contribution from the x^6 term    = {C0_term2:.6f}")
    print(f"Contribution from the x^8 term    = {C0_term3:.6f}")


solve_integral_approximation()
<<<I(epsilon) approx = 0.6888 * epsilon^(-4/5) - 0.0754 * epsilon^(-3/5) - 0.1488 * epsilon^(-1/5) + 0.0074>>>