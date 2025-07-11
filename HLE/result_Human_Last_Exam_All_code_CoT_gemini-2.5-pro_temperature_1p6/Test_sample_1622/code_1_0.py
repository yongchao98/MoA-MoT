def solve():
    """
    This function prints the derived formula for P(n).
    """
    # The derived formula for P(n) where L = ln(n)
    # P(n) is composed of terms to correct the approximation up to the required error bound.
    p_n_formula = "( (3*L**2 + 2*L - 2) / (24*n**2) ) + ( (L**3 + 2*L**2 - 2*L) / (48*n**3) )"

    print("The formula for P(n) is:")
    print("P(n) = " + p_n_formula)

solve()