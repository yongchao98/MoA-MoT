import sympy

def solve_and_print_expression():
    """
    This function derives and prints the exact value of l_k(n).
    The derivation is based on the properties of the Hyperbolic Normal distribution.
    The final expression is formulated and printed using the sympy library.
    """
    
    # Declare the variables n and k as symbolic variables
    n, k = sympy.symbols('n k', real=True)

    # Based on the derivation, we define the three main components of the expression for l_k(n)
    
    # Term 1: From the determinant of the covariance matrix
    # det(Sigma) = 1 / (n + 1)
    # -1/2 * log(det(Sigma)) = 1/2 * log(n + 1)
    term1_val = sympy.Rational(1, 2)
    term1_log = sympy.log(n + 1)
    term1 = term1_val * term1_log
    
    # Term 2: From the quadratic form in the exponent of the Gaussian density
    # -1/2 * n^T * Sigma^-1 * n = -1/2 * (k**2/n) * (4*n - 2) = -k**2 * (2 - 1/n)
    term2_val_n = sympy.Rational(1, n)
    term2 = -k**2 * (2 - term2_val_n)
    
    # Term 3: From the Jacobian of the exponential map
    # (n - 1) * log(k/sinh(k)), where sinh(k) = 1
    # (n - 1) * log(k)
    term3_log = sympy.log(k)
    term3 = (n - 1) * term3_log
    
    # The complete expression for l_k(n)
    l_k_n_expression = term1 + term2 + term3
    
    # Print the final result in a clear, equation-like format
    print("The exact value of l_k(n) is:")
    sympy.pprint(sympy.Eq(sympy.Symbol('l_k(n)'), l_k_n_expression, evaluate=False), use_unicode=True)

solve_and_print_expression()