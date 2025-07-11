import sympy as sp

def solve_lk_n():
    """
    Calculates and prints the exact symbolic expression for l_k(n).
    """
    # Define symbolic variables for n and k.
    # We assume n is an integer >= 3, and k is a constant.
    n, k = sp.symbols('n k', real=True)

    # Derived expression for l_k(n) based on the step-by-step derivation.
    # The expression consists of three parts:
    # 1. The term from the determinant of Sigma.
    term_det = sp.log(n + 1) / 2
    # 2. The term from the quadratic form n^T * Sigma^-1 * n.
    term_quad = -k**2 * (2*n - 1) / n
    # 3. The term from the volume element correction factor of the hyperbolic space.
    term_jacobian = (n - 1) * sp.log(k)

    # Combine the terms to get the final expression for l_k(n).
    l_k_n_expression = term_det + term_quad + term_jacobian

    # Print the final result in a clear format.
    print("The exact value of l_k(n) is:")
    
    # We construct a formatted string to show the final equation clearly,
    # as requested by the prompt ("output each number in the final equation").
    expression_string = f"1/2*log(n + 1) - k**2*(2*n - 1)/n + (n - 1)*log(k)"
    
    # You can also print the sympy object directly:
    # print(l_k_n_expression)

    print(expression_string)

solve_lk_n()