def solve():
    """
    This function prints the derived exact formula for l_k(n).
    The derivation is based on the properties of the Hyperbolic Normal distribution
    and the specific parameters mu and Sigma provided in the problem.
    """

    # The final derived formula for l_k(n) is composed of three terms.
    # We construct the string representing this formula.
    # The numbers in the equation (1, 2, 1, 1) are explicitly included.
    
    formula = "l_k(n) = (1/2)*ln(n + 1) - k^2*((n + 1)/n) + (n - 1)*ln(k)"
    
    print(formula)

solve()