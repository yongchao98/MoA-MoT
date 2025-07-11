import math

def calculate_performance_bound(H, abs_A, lambda_val):
    """
    Calculates the tightest upper bound for J(pi^*) - J(hat_pi).

    The performance difference between the expert policy and the learned policy
    in imitation learning can be bounded. A standard result for a fixed learned
    policy (as in Behavioral Cloning) shows that the difference accumulates
    quadratically with the horizon H.

    The bound is derived from the performance difference lemma:
    J(pi^*) - J(hat_pi) <= sum_{t=0 to H-1} E_{s_t~d_pi_hat} [ (H-t) * TV(pi*(.|s_t), hat_pi(.|s_t)) ]

    Using the given upper bound on the total variation risk:
    TV(pi*, hat_pi) <= |A| * (1 - e^(-lambda))

    The final upper bound is: H * (H+1) / 2 * |A| * (1 - e^(-lambda))

    Args:
        H (int): The episode horizon.
        abs_A (int): The size of the discrete action space |A|.
        lambda_val (float): The hyperparameter lambda.

    Returns:
        float: The tightest upper bound on the performance difference.
    """
    if H <= 0 or abs_A <= 0 or lambda_val < 0:
        raise ValueError("Parameters H, abs_A must be positive, and lambda_val must be non-negative.")

    tv_risk_bound = abs_A * (1 - math.exp(-lambda_val))
    
    # The sum of (H-t) from t=0 to H-1 is H + (H-1) + ... + 1 = H*(H+1)/2
    horizon_factor = H * (H + 1) / 2
    
    upper_bound = horizon_factor * tv_risk_bound
    
    # Printing the breakdown of the calculation as requested.
    # The reward range R_max is assumed to be 1 as it's not specified.
    print(f"The performance difference J(pi^*) - J(hat_pi) is upper bounded by:")
    print(f"Bound = [H * (H + 1) / 2] * |A| * (1 - exp(-lambda))")
    print(f"      = [{H} * ({H} + 1) / 2] * {abs_A} * (1 - exp(-{lambda_val:.4f}))")
    print(f"      = [{horizon_factor}] * {abs_A} * ({1 - math.exp(-lambda_val):.4f})")
    print(f"      = {upper_bound}")

    return upper_bound

# Example usage with some plausible parameters.
# H: Horizon of the task.
# abs_A: Number of possible actions.
# lambda_val: A hyperparameter from the learning algorithm.
H_example = 20
abs_A_example = 5
lambda_example = 0.05

calculate_performance_bound(H_example, abs_A_example, lambda_example)