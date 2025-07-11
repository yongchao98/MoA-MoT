def solve_hamiltonicity_threshold():
    """
    This function constructs and prints the symbolic formula for the d-threshold of Hamiltonicity.
    """

    # The problem asks for the d-threshold for Hamiltonicity for H_n U G(n, p),
    # where the minimum degree d is given by d = n/2 - eta.

    # Based on known results for the related problem of a random graph G(n,p)
    # conditioned on minimum degree, a precise formula for the threshold p is:
    # p = (ln(n/d) + ln(ln(n/d))) / d
    # We will use this formula and substitute d = n/2 - eta.

    # Define the variables as strings for symbolic representation
    n_str = "n"
    eta_str = "η"

    # The only number in the formula is 2.
    number_two = 2

    # Create the expression for d
    d_expression = f"({n_str}/{number_two} - {eta_str})"

    # Create the expression for the term n/d
    n_over_d_expression = f"{n_str} / {d_expression}"

    # Build the numerator of the final formula
    numerator_expression = f"ln({n_over_d_expression}) + ln(ln({n_over_d_expression}))"

    # Build the denominator of the final formula
    denominator_expression = d_expression

    # Construct the full equation as a string
    final_equation = f"p = ({numerator_expression}) / ({denominator_expression})"

    print("The d-threshold for Hamiltonicity in the given range is given by the formula:")
    print(final_equation)

solve_hamiltonicity_threshold()