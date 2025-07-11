def solve_asymptotic_degree():
    """
    This function prints the asymptotic formula for d_{B,delta}.
    The formula is derived by analyzing the behavior in different regimes of B and L,
    where L = log(1/delta).

    The final expression combines the behavior from the "large interval" regime (B)
    and the "high precision" regime (L/log(L)).
    """

    # The formula is B + L / log(L).
    # Following the instruction to output each number, we format the string this way.
    # The coefficients and constants are all implicitly 1.
    final_equation = "1*B + 1*L / log(1*L)"
    
    print(final_equation)

solve_asymptotic_degree()