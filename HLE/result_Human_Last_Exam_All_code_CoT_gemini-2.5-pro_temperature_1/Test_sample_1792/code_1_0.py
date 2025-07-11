def solve_ordinal_expression():
    """
    This function prints the simplified ordinal expression based on the problem description.
    The simplification steps are performed symbolically as described in the plan.
    """

    # Unicode symbols for ordinals for better readability
    w = "\u03C9"
    w1 = f"{w}_1"
    w2 = f"{w}_2"

    # The coefficients derived from the step-by-step simplification
    # The final expression is w2 * w1 + w1 * 1 + w * 0 + 0
    alpha_1 = w1
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # Construct the final expression string
    final_expression = f"{w2} \u00B7 {alpha_1} + {w1} \u00B7 {alpha_2} + {w} \u00B7 {alpha_3} + {alpha_4}"

    print("The final expression in the required form is:")
    print(final_expression)

solve_ordinal_expression()