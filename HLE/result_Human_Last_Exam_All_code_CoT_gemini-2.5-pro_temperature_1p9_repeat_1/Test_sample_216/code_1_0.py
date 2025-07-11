def solve_imitation_learning_bound():
    """
    This function formulates and prints the tightest upper bound for the performance
    difference J(pi^*) - J(pi_hat) based on the provided information.
    """

    # Define the symbolic components of the equation
    performance_difference = "J(pi^*) - J(pi_hat)"
    horizon = "H"
    max_reward = "R_max" # Maximum single-step reward
    action_space_size = "|A|"
    hyperparameter = "lambda"

    # The equation for the upper bound is derived from standard imitation learning theory.
    # The performance difference is bounded by H^2 * R_max * (population_TV_risk).
    # We substitute the given upper bound for the population_TV_risk.

    # Print the components and the final inequality
    print("Based on the standard compounding error analysis in imitation learning:")
    print(f"The performance difference is bounded as follows:")
    
    bound_expression = f"{horizon}^2 * {max_reward} * {action_space_size} * (1 - exp(-{hyperparameter}))"
    
    # Print the final inequality with all terms
    print(f"{performance_difference} <= {bound_expression}")
    
    print("\nWhere:")
    print(f"  {performance_difference}: The difference in expected return between the expert and the learned policy.")
    print(f"  {horizon}: The episode horizon.")
    print(f"  {max_reward}: The maximum possible single-step reward.")
    print(f"  {action_space_size}: The size of the action space.")
    print(f"  {hyperparameter}: A hyperparameter from the learning algorithm.")
    print(f"  exp(): The exponential function.")
    
solve_imitation_learning_bound()