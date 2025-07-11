import numpy as np
from scipy.integrate import quad

def analyze_functions():
    """
    Analyzes the inequality for the three given functions.
    """
    print("Analyzing the inequality Sum(n*|a_n|^2) <= Sum(|a_n|) for each function.\n")

    # --- Function 1 ---
    print("--- Function 1: f(z) = sum(z^(2^(2^n)) / 2^n) ---")
    rhs_1 = 0
    for n in range(1, 50): # The sum to infinity is 1
        rhs_1 += 1/2**n
    
    lhs_terms_1 = []
    for n in range(1, 5):
        term = (2**(2**n)) / (2**(2*n))
        lhs_terms_1.append(term)
    
    print(f"RHS = sum(1/2^n) = {rhs_1}")
    print(f"LHS = sum( (2^(2^n)) / (2^(2n)) )")
    print(f"The first few terms of the LHS sum are: {lhs_terms_1[0]}, {lhs_terms_1[1]}, {lhs_terms_1[2]}, {lhs_terms_1[3]}, ...")
    total_sum_4_terms = sum(lhs_terms_1)
    print(f"The sum of just the first 4 terms is {total_sum_4_terms}. The series diverges to infinity.")
    print("Equation: infinity <= 1")
    print("Conclusion for 1: The inequality is FALSE.\n")

    # --- Function 2 ---
    print("--- Function 2: Conformal map to a square ---")
    # S = 2 * Integral from 0 to 1 of (1-x^4)^(-1/2) dx
    integrand_S = lambda x: (1 - x**4)**(-0.5)
    S, _ = quad(integrand_S, 0, 1)
    S *= 2
    
    # S2 = 2 * Integral from 0 to 1 of (1+x^4)^(-1/2) dx
    integrand_S2 = lambda x: (1 + x**4)**(-0.5)
    S2, _ = quad(integrand_S2, 0, 1)
    S2 *= 2
    
    # LHS = Area/pi = |C|^2 * S^2 / pi, where |C|^2 = 2
    lhs_2 = 2 * S**2 / np.pi
    
    # RHS = |a_0| + |C| * sum(|g_n|) = S2 + sqrt(2) * (S/2)
    rhs_2 = S2 + np.sqrt(2) * S / 2
    
    print(f"LHS = 2 * ({S:.4f})^2 / pi = {lhs_2:.4f}")
    print(f"RHS = {S2:.4f} + sqrt(2) * {S:.4f} / 2 = {rhs_2:.4f}")
    print(f"Equation: {lhs_2:.4f} <= {rhs_2:.4f}")
    print("Conclusion for 2: The inequality is FALSE.\n")

    # --- Function 3 ---
    print("--- Function 3: Conformal map to the Koch snowflake interior ---")
    print("LHS = Area(Koch snowflake) / pi. The area is finite, so LHS is finite.")
    print("RHS = sum(|a_n|). This sum diverges to infinity because the boundary of the Koch snowflake is not rectifiable (has infinite length).")
    print("Equation: [a finite number] <= infinity")
    print("Conclusion for 3: The inequality is TRUE.\n")

analyze_functions()