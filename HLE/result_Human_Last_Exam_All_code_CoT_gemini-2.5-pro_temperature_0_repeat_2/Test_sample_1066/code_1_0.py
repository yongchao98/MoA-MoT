def print_lower_bound_formula():
    """
    This function prints the derived lower bound for E[S].
    The formula involves the number of tokens n, the average entropy α,
    and the constants 1, 6, and π.
    """
    n_symbol = "n"
    alpha_symbol = "α"
    pi_symbol = "π"
    
    number_one = 1
    number_six = 6
    number_two = 2

    # Construct and print the formula string
    formula = f"{n_symbol} * ({alpha_symbol} + {number_one} - {number_six}/({pi_symbol}^{number_two}))"
    print("A lower bound on E[S] is given by the expression:")
    print(formula)

print_lower_bound_formula()