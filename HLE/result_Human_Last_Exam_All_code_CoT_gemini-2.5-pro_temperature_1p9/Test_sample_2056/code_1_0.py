def solve():
    """
    This function prints the derived exact value of l_k(n).
    The derivation involves evaluating the log-probability density function of a
    wrapped normal distribution on the hyperboloid at a specific point.
    The final expression is found by algebraically simplifying the components
    of the PDF formula, which include the hyperbolic distance, the inverse
    sampling map, and the determinant and inverse of the covariance matrix.
    """

    # The final derived expression for l_k(n) has three parts:
    # 1. A term from the determinant of Sigma: (1/2) * ln(n+1)
    # 2. A term from the Jacobian of the exponential map: (n-1) * ln(k)
    # 3. A term from the Mahalanobis distance in the normal PDF exponent: -k^2 * (2 - 1/n)
    
    # We construct the string to represent the final formula clearly.
    # The instructions require printing each number in the equation.
    term1_coeff = "1/2"
    term1_arg = "n + 1"
    
    term2_coeff = "n - 1"
    
    term3_factor1_exp = "2"
    term3_factor2_term1 = "2"
    term3_factor2_term2 = "1/n"
    
    final_expression = (f"({term1_coeff}) * ln({term1_arg}) + "
                        f"({term2_coeff}) * ln(k) - "
                        f"k**{term3_factor1_exp} * ({term3_factor2_term1} - {term3_factor2_term2})")
    
    print(f"l_k(n) = {final_expression}")

solve()