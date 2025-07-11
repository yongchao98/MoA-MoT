def solve():
    """
    This function prints the two additional inequalities needed to complete the MILP model.
    """
    # The first inequality is y >= x - a - b
    # The coefficients for x, a, b are 1, -1, and -1 respectively.
    inequality1 = "y >= x - a - b"

    # The second inequality is y >= l - l*b
    # The parameter 'l' is the lower bound for x.
    inequality2 = "y >= l - l*b"

    print(inequality1 + ", " + inequality2)

solve()