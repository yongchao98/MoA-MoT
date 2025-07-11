import math

def solve_imitation_learning_bound():
    """
    This function calculates and prints the symbolic formula for the tightest upper bound
    of the performance difference J(pi^*) - J(pi_hat) based on the provided information.

    The derivation relies on standard imitation learning theory, which results in a bound
    that depends quadratically on the horizon H.
    """

    # The formula for the upper bound is derived as H^2 * |A| * (1 - e^(-lambda)).
    # The components of this formula are:
    # H: The horizon of the episode.
    # |A|: The size of the discrete action space.
    # lambda: A hyperparameter of the algorithm.

    # According to the prompt, we need to output each number in the final equation.
    # The numbers in this formula are:
    # 2: The exponent for the horizon H.
    # 1: The constant term inside the parenthesis.
    # -1: The coefficient of lambda in the exponent of e.

    power_of_H = 2
    constant_term = 1
    lambda_coefficient = -1

    # We use string representations for the symbolic parts of the formula.
    H_symbol = "H"
    A_size_symbol = "|A|"
    lambda_symbol = "lambda"

    # Construct the final formula string.
    final_formula = (
        f"{H_symbol}^{power_of_H} * {A_size_symbol} * "
        f"({constant_term} - exp({lambda_coefficient} * {lambda_symbol}))"
    )

    print("The tightest upper bound of J(pi^*) - J(pi_hat) is:")
    print(final_formula)

solve_imitation_learning_bound()