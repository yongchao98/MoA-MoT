def solve():
    """
    This function determines and prints the expression for the upper bound H.
    """
    # The parameters are symbolic. We will construct the formula as a string.
    # a = k, b = ||rho(0,.)||_L1, c = pi, d = nu, r = lower bound of rho, t = time
    # The derived upper bound H is (-a * b * t) / (c * d^2 * r).
    # We explicitly write out the coefficients as requested by the prompt.

    H_expression = "H(a, b, c, d, r, t) = (-1 * a * b * t) / (c * d**2 * r)"
    print(H_expression)

solve()