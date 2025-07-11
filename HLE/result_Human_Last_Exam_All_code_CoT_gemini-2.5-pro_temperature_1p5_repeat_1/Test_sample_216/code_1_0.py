import math

def calculate_performance_gap_bound(H, A_size, lambda_val):
    """
    Calculates and prints the tightest upper bound on the imitation learning performance gap
    J(pi^*) - J(pi_hat) based on the provided TV risk.

    Args:
        H (int): The episode horizon.
        A_size (int): The size of the action space, |A|.
        lambda_val (float): The hyperparameter lambda.
    """
    # The performance gap J(pi*) - J(pi_hat) is bounded by H^2 * T(pi_hat, pi*).
    # We are given T(pi_hat, pi*) <= |A| * (1 - e^(-lambda)).
    # Therefore, the tightest upper bound is H^2 * |A| * (1 - e^(-lambda)).

    # Calculate the three numerical components of the final equation
    horizon_factor = float(H**2)
    action_space_size = float(A_size)
    risk_factor = 1.0 - math.exp(-lambda_val)

    # Calculate the final upper bound
    upper_bound = horizon_factor * action_space_size * risk_factor

    # As requested, output each number in the final equation, then the full equation and result.
    print("The final equation for the upper bound is of the form: (H^2) * |A| * (1 - exp(-lambda))")
    print("\nThe numerical values of the components are:")
    print(f"H^2 = {horizon_factor}")
    print(f"|A| = {action_space_size}")
    print(f"1 - exp(-lambda) = {risk_factor}")

    print("\n---")
    print("Final equation with values:")
    # Using 'g' format to avoid excessive decimals for the risk factor
    print(f"{horizon_factor} * {action_space_size} * {risk_factor:.5g} = {upper_bound}")


# --- Example Usage ---
# You can change these values to match your specific problem.
# H: The episode horizon. Must be >= 0.
H_example = 50

# A_size: The number of discrete actions |A|. Must be > 0.
A_size_example = 4

# lambda_val: A positive hyperparameter of the algorithm.
lambda_example = 0.1

calculate_performance_gap_bound(H_example, A_size_example, lambda_example)