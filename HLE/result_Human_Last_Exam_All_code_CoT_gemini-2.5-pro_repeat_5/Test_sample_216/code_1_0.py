def solve_imitation_learning_bound():
    """
    This function formulates and prints the tightest upper bound for J(pi^*) - J(pi_hat).
    The derivation leads to the formula: 2 * H * |A| * (1 - exp(-lambda)).
    """
    
    # Define variables as string representations for the formula
    H = "H"
    abs_A = "|A|"
    lambda_val = "lambda"
    
    # The final expression for the tightest upper bound
    # It contains the number 2, and implicitly 1 and -1 inside the exponential term.
    upper_bound_expression = f"2 * {H} * {abs_A} * (1 - exp(-{lambda_val}))"
    
    print("The tightest upper bound for the performance difference J(pi^*) - J(pi_hat) is:")
    print(upper_bound_expression)

solve_imitation_learning_bound()