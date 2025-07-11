import math

def calculate_performance_difference_bound(H, A_size, lambda_val):
    """
    Calculates the tightest upper bound of J(pi*) - J(pi_hat).

    Args:
        H (int): The horizon of the episode.
        A_size (int): The size of the discrete action space |A|.
        lambda_val (float): The hyperparameter lambda.
    """
    # The population total variation (TV) risk is bounded by:
    # TV_risk <= |A| * (1 - e^(-lambda))
    tv_risk_bound = A_size * (1 - math.exp(-lambda_val))

    # The performance difference J(pi*) - J(pi_hat) is bounded by:
    # J(pi*) - J(pi_hat) <= H^2 * TV_risk
    performance_diff_bound = H**2 * tv_risk_bound

    # Printing the results and the equation
    print("Given the parameters:")
    print(f"  Horizon (H) = {H}")
    print(f"  Action space size (|A|) = {A_size}")
    print(f"  Hyperparameter (lambda) = {lambda_val}\n")

    print("The performance difference J(pi*) - J(pi_hat) has a tight upper bound given by the equation:")
    print(f"  J(pi*) - J(pi_hat) <= H^2 * |A| * (1 - e^(-lambda))")
    print("Plugging in the values:")
    print(f"  J(pi*) - J(pi_hat) <= {H}^2 * {A_size} * (1 - e^(-{lambda_val}))")
    print(f"  J(pi*) - J(pi_hat) <= {H**2} * {A_size} * {1 - math.exp(-lambda_val):.4f}")
    print(f"  J(pi*) - J(pi_hat) <= {performance_diff_bound:.4f}\n")
    
    print("Final computed bound:")
    print(performance_diff_bound)

if __name__ == '__main__':
    # Example values for the parameters
    H = 20         # Episode horizon
    A_size = 5     # Size of the action space
    lambda_val = 0.1 # Hyperparameter
    
    calculate_performance_difference_bound(H, A_size, lambda_val)
