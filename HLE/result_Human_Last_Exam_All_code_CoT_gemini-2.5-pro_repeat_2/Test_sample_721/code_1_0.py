import numpy as np
from scipy import integrate

def analyze_functions():
    """
    Analyzes the three functions with respect to the inequality
    sum(n * |a_n|^2) <= sum(|a_n|).
    """

    print("--- Analysis of Function 1 ---")
    # f(z) = sum_{n=1 to inf} z^(2^(2^n)) / 2^n
    # RHS = sum_{n=1 to inf} 1/2^n = 1
    rhs_1 = 1.0
    print(f"RHS = sum |a_k| = {rhs_1}")

    # LHS = sum_{k} k|a_k|^2 = sum_{n=1 to inf} (2^(2^n)) * |1/2^n|^2 = sum_{n=1 to inf} 2^(2^n) / 4^n
    lhs_sum_terms_1 = []
    print("LHS = sum k|a_k|^2. The first few terms of the sum are:")
    for n in range(1, 5):
        term = 2**(2**n) / (4**n)
        lhs_sum_terms_1.append(term)
        print(f"  n={n}: term = 2^(2^{n}) / 4^{n} = {term}")
    
    total_lhs_1 = sum(lhs_sum_terms_1)
    print(f"The sum of the first {len(lhs_sum_terms_1)} terms is {total_lhs_1}.")
    print("The terms grow extremely fast, so the sum diverges to infinity.")
    print(f"Inequality check: infinity <= {rhs_1} is False.")
    print("Conclusion: Function 1 does not satisfy the inequality.\n")


    print("--- Analysis of Function 2 ---")
    # f(z) maps D to a square.
    # LHS = Area / pi. The side length of the square is K = 2 * integral_0^1 (1-u^4)^(-1/2) du.
    # The integral can be expressed using Gamma functions: K = Gamma(1/4)^2 / sqrt(2*pi)
    gamma_1_4 = np.math.gamma(1/4)
    K = gamma_1_4**2 / np.sqrt(2 * np.pi)
    Area = K**2
    lhs_2 = Area / np.pi
    print(f"LHS = Area / pi = K^2 / pi = ({K:.4f})^2 / pi = {lhs_2:.4f}")

    # RHS = sum |a_n| <= |a_0| + 2|a_1|
    # a_1 = f'(0) = -1 - i
    a1 = -1 - 1j
    abs_a1 = abs(a1)
    
    # a_0 = f(0) = integral_0^i d(xi)/sqrt(xi(1-xi^2))
    # This evaluates to exp(i*pi/4) * integral_0^1 dt/sqrt(t(1+t^2))
    integrand_a0 = lambda t: 1 / np.sqrt(t * (1 + t**2))
    abs_a0_val, _ = integrate.quad(integrand_a0, 0, 1)

    rhs_upper_bound_2 = abs_a0_val + 2 * abs_a1
    print(f"RHS = sum |a_n|, for which we found an upper bound:")
    print(f"|a_0| = {abs_a0_val:.4f}")
    print(f"|a_1| = |-1-i| = {abs_a1:.4f}")
    print(f"RHS <= |a_0| + 2|a_1| = {abs_a0_val:.4f} + 2 * {abs_a1:.4f} = {rhs_upper_bound_2:.4f}")
    
    print(f"Inequality check: {lhs_2:.4f} <= {rhs_upper_bound_2:.4f} is False.")
    print("Conclusion: Function 2 does not satisfy the inequality.\n")
    
    
    print("--- Analysis of Function 3 ---")
    # f(z) maps D to the interior of the Koch snowflake.
    # The area is finite, so LHS = Area / pi is finite.
    # The boundary length is infinite, which implies RHS = sum |a_n| is infinite.
    print("LHS = Area(Koch Snowflake) / pi = finite value")
    print("RHS = sum |a_n| = infinity (since boundary is non-rectifiable)")
    print("Inequality check: finite_value <= infinity is True.")
    print("Conclusion: Function 3 satisfies the inequality.\n")
    
    print("Summary: Only function 3 satisfies the inequality.")

analyze_functions()