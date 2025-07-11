def solve_work_cycle():
    """
    This function derives and prints the formula for the work done by the current source.
    The derivation steps are outlined in the explanation. The code constructs the final
    formula string based on the derived components.
    """

    # Step 1: Define the components of the formula based on the derivation.
    # The change in inductance, L(x2) - L(x1), simplifies to:
    # (N^2 * w / g) * (mu - mu_0) * (x2 - x1)
    
    # The total work W simplifies to:
    # W = -1/2 * (L(x2) - L(x1)) * (I2^2 - I1^2)

    # We will build the final expression from its parts for clarity.
    
    # Numerator part of the main fraction
    numerator_terms = [
        "(mu - mu_0)",
        "N^2",
        "w",
        "(x2 - x1)",
        "(I2^2 - I1^2)"
    ]
    
    # Denominator part of the main fraction
    denominator_terms = ["2", "g"]

    # Combine the terms to form the final expression string
    # We lead with a negative sign from the derivation W = -1/2 * ...
    final_equation = "W = - "
    
    # Assemble the fraction
    # In many physics notations, the constants are grouped into a single fraction
    final_equation += f"({numerator_terms[0]}) / ({' * '.join(denominator_terms)})"
    
    # Append the remaining terms
    final_equation += f" * {numerator_terms[1]}"
    final_equation += f" * {numerator_terms[2]}"
    final_equation += f" * {numerator_terms[3]}"
    final_equation += f" * {numerator_terms[4]}"

    print("The work done by the current source for each cycle is calculated as follows:")
    
    # Print the final equation with each variable and number explicitly shown.
    # W = - (mu - mu_0) / (2*g) * N^2 * w * (x2 - x1) * (I2^2 - I1^2)
    # The following print statement produces a format consistent with the derivation result
    # and the multiple-choice options.
    
    print("W = - (mu - mu_0) * N^2 * w * (x2 - x1) * (I2^2 - I1^2) / (2 * g)")

    # For better readability, matching the structure of the options:
    print("\nOr, written in a form similar to the choices:")
    
    mu_diff_term = "(mu - mu_0)"
    factor_2g_inv = f"(2*g)"
    N_sq_term = "N^2"
    w_term = "w"
    x_diff_term = "(x2 - x1)"
    I_sq_diff_term = "(I2^2 - I1^2)"

    final_formatted_eq = f"W = - ({mu_diff_term} / {factor_2g_inv}) * {N_sq_term} * {w_term} * {x_diff_term} * {I_sq_diff_term}"
    
    # To avoid messy division in the first term for the final output, let's rearrange to match answer D
    final_formatted_eq = f"W = - (({mu_diff_term}) / (2*g)) * {N_sq_term} * {w_term} * {x_diff_term} * {I_sq_diff_term}"

    # Re-arranging for final clear display matching option D.
    final_equation_string = f"W = - ( (mu - mu_0) / (2 * g) ) * N^2 * w * (x2 - x1) * (I2^2 - I1^2)"
    
    # To output each component as requested:
    term1 = "-(mu - mu_0)"
    term2 = "2*g"
    term3 = "N^2"
    term4 = "w"
    term5 = "(x2 - x1)"
    term6 = "(I2^2 - I1^2)"

    print(f"Final equation: W = ( {term1} / ({term2}) ) * {term3} * {term4} * {term5} * {term6}")

solve_work_cycle()