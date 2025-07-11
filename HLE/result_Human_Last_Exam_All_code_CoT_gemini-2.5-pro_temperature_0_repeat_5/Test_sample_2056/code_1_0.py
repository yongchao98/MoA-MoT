def print_exact_formula():
    """
    This function prints the final derived formula for l_k(n).
    The derivation involves analyzing the hyperbolic normal distribution defined
    by the sampling procedure and evaluating its log-pdf at the specified point.
    The final expression is given in terms of n and k.
    """

    # The derived formula is composed of three parts:
    # 1. From the determinant of Sigma: (n + 1)/2 * ln(n + 1)
    # 2. From the Jacobian of the exponential map: (n - 1) * ln(k)
    # 3. From the quadratic form in the exponent: -k**2 * (n + 1) * (2*n - 1) / n
    
    # The problem asks to output the numbers in the final equation.
    # The following string represents the equation with all numbers explicitly written.
    
    final_formula = "(n + 1)/2 * ln(n + 1) + (n - 1) * ln(k) - k**2 * (n + 1) * (2*n - 1) / n"
    
    print("The exact value of l_k(n) is given by the following expression:")
    print(final_formula)

print_exact_formula()