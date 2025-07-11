import sympy as sp

def solve_lk_n_formula():
    """
    This function derives and prints the symbolic formula for l_k(n).
    The derivation is based on the analysis of the probability density function
    of the given sampling procedure evaluated at a specific point.
    """
    
    # Define symbols for n and k, assuming n>=3 is an integer, k is a real constant.
    n, k = sp.symbols('n k', real=True, positive=True)

    # The formula for l_k(n) is derived as a sum of several terms.
    # l_k(n) = ln[p(x)] + (n/2)*ln(2*pi)
    # where ln[p(x)] = - (n/2)*ln(2*pi) - (1/2)*ln(det(Sigma)) - ln(mu_1)
    #                   + (n-1)*ln(d/sinh(d)) - (1/2)*n^T*Sigma^-1*n
    #
    # After calculation, the terms are:
    # det(Sigma) = 1/(n+1)
    # mu_1 = sqrt(2)
    # d = k, sinh(d) = 1
    # n^T*Sigma^-1*n = k^2 * (4 - 2/n)
    #
    # So, l_k(n) = -(1/2)*ln(1/(n+1)) - ln(sqrt(2)) + (n-1)*ln(k) - (1/2)*k^2*(4-2/n)
    # Which simplifies to the expression below.

    # Term 1: From the determinant of the covariance matrix Sigma
    term1 = sp.Rational(1, 2) * sp.log(n + 1)
    
    # Term 2: From the factor 1/mu_1 in the PDF
    term2 = -sp.Rational(1, 2) * sp.log(2)

    # Term 3: From the Jacobian of the exponential map
    term3 = (n - 1) * sp.log(k)
    
    # Term 4: From the Mahalanobis distance in the exponent of the normal distribution
    term4 = -(2 - sp.Rational(1, n)) * k**2

    # The complete symbolic expression for l_k(n)
    l_k_n_expression = term1 + term2 + term3 + term4

    # We print the final formula, showing each number and component term by term,
    # as requested in the instructions.
    
    print("The exact formula for l_k(n) is given by the sum of the following terms:")
    print(f"Term from det(Sigma): {term1}")
    print(f"Term from mean's first component: {term2}")
    print(f"Term from the Jacobian: {term3}")
    print(f"Term from the Gaussian exponent: {term4}")
    
    print("\nThus, the final expression for l_k(n) is:")
    # Using sp.pretty_print for a clear mathematical representation of the final formula.
    sp.pprint(l_k_n_expression, use_unicode=True)

if __name__ == '__main__':
    solve_lk_n_formula()
