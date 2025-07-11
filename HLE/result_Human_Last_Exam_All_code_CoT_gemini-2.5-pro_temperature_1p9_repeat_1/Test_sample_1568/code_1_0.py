def solve_and_print_product_formula():
    """
    This function generates and prints the symbolic result for the infinite product:
    Product_{n=3 to infinity} (1 - z^3 / n^3)

    The code constructs the final equation from its numerical components as requested.
    """

    # Define the numbers present in the problem statement and the final equation.
    n_start = 3
    power = 3
    term1_n = 1
    term2_n = 2
    term2_denominator = term2_n**power

    # Build the Left Hand Side (LHS) of the equation.
    lhs_product_symbol = "\u220F"
    infinity_symbol = "\u221E"
    lhs = f"{lhs_product_symbol}_{{n={n_start}}}^{infinity_symbol}(1 - z^{power}/n^{power})"

    # Build the Right Hand Side (RHS) of the equation.
    numerator = term2_denominator

    # Denominator components
    denominator_factor1 = f"(1 - z^{power})"
    denominator_factor2 = f"({term2_denominator} - z^{power})"
    
    # Gamma function components
    gamma_symbol = "\u0393"
    pi_symbol = "\u03C0"
    
    # Arguments of the Gamma functions
    gamma_arg1 = "1 - z"
    gamma_arg2 = f"1 - z*exp(2{pi_symbol}i/{power})"
    gamma_arg3 = f"1 - z*exp(4{pi_symbol}i/{power})"

    gamma_part = f"{gamma_symbol}({gamma_arg1}){gamma_symbol}({gamma_arg2}){gamma_symbol}({gamma_arg3})"
    
    rhs = f"{numerator} / ({denominator_factor1}{denominator_factor2}{gamma_part})"

    # Print the final, fully assembled equation.
    print("The final expression for the infinite product is:")
    print(f"{lhs} = {rhs}")

solve_and_print_product_formula()