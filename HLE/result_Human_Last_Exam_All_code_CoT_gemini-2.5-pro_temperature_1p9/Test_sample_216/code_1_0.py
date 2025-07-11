import math

def calculate_imitation_learning_bound(H_symbol, A_size_symbol, lambda_symbol):
    """
    Calculates and prints the symbolic upper bound for the performance difference
    in an imitation learning algorithm.

    Args:
        H_symbol (str): The symbol for the episode horizon.
        A_size_symbol (str): The symbol for the size of the action space.
        lambda_symbol (str): The symbol for the hyperparameter lambda.
    """
    # The derivation leads to the bound: 2 * H^2 * Risk
    # We are given Risk <= |A| * (1 - exp(-lambda))
    # So the final bound is 2 * H^2 * |A| * (1 - exp(-lambda))

    # Constants and parts of the equation
    constant_factor = "2"
    horizon_term = f"{H_symbol}^2"
    action_space_term = f"{A_size_symbol}"
    lambda_term = f"(1 - e^(-{lambda_symbol}))"

    # Assemble the final expression
    final_bound_expression = f"{constant_factor} * {horizon_term} * {action_space_term} * {lambda_term}"

    # Print the result clearly, showing each number and symbol
    print("The tightest upper bound of J(pi^*) - J(pi_hat) is:")
    print(final_bound_expression)
    
    # You can uncomment the following lines to compute a numerical value for specific inputs
    # H = 10
    # A_size = 4
    # lmbda = 0.1
    # numeric_bound = 2 * (H**2) * A_size * (1 - math.exp(-lmbda))
    # print(f"\nFor example, if H={H}, |A|={A_size}, and lambda={lmbda}:")
    # print(f"Bound = {numeric_bound}")

if __name__ == "__main__":
    # Define the symbols as strings for the final expression
    H_symbol = "H"
    A_size_symbol = "|A|"
    lambda_symbol = "lambda"
    
    calculate_imitation_learning_bound(H_symbol, A_size_symbol, lambda_symbol)
