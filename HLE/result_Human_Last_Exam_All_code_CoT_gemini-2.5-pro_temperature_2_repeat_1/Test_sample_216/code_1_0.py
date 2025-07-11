import math

def display_performance_bound():
    """
    This function formulates and prints the tightest upper bound for the
    performance difference J(pi^*) - J(pi_hat) based on the provided
    Total Variation (TV) risk in an imitation learning setting.

    The analysis relies on standard results from imitation learning theory,
    specifically the quadratic dependence of the performance gap on the horizon H
    in the behavioral cloning setup due to compounding errors.
    """
    # Let H be the episode horizon.
    # Let |A| be the size of the action space.
    # Let lambda be the hyperparameter.

    # A standard theoretical bound on the performance difference is:
    # J(pi^*) - J(pi_hat) <= (H - 1)^2 * R_max * epsilon
    # where epsilon is the one-step TV risk max_s d_TV(pi^*, pi_hat).
    # Assuming R_max = 1 as it's not specified.

    # The problem provides the bound on the TV risk:
    # T(pi_hat, pi_star) <= |A| * (1 - exp(-lambda))
    
    # We substitute this into the performance difference bound to get the final result.
    
    # The question requests that all numbers in the final equation be output.
    # The numbers are: 1 (in H-1), 2 (the exponent), 1 (in 1-exp), and -1 (in -lambda).
    final_equation = "J(pi^*) - J(pi_hat) <= (H - 1)^2 * |A| * (1 - exp(-1 * lambda))"
    
    print("Based on standard imitation learning theory for the Behavioral Cloning setting,")
    print("the tightest upper bound on the performance difference is given by the following equation:")
    print(final_equation)

# Execute the function to display the result.
display_performance_bound()