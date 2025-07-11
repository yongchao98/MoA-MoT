import math

def solve_imitation_learning_bound():
    """
    This function calculates and prints the symbolic expression for the tightest
    upper bound of the performance difference J(pi^*) - J(pi_hat).
    """
    # The problem provides symbolic variables, not numerical values.
    # We represent them as strings for the final expression.
    H = "H"  # Episode horizon
    A_size = "|A|"  # Size of the discrete action space
    lam = "lambda"  # Hyperparameter from the algorithm

    # The performance difference J(pi^*) - J(pi_hat) is bounded by H^2 * TV_risk.
    # This is a standard result from imitation learning theory, where the error
    # compounds quadratically with the horizon H.

    # We are given the upper bound for the TV risk:
    # TV_risk <= |A| * (1 - exp(-lambda))

    # To find the tightest upper bound on the performance difference, we substitute
    # the upper bound of the TV risk into the performance gap formula.
    # The resulting expression is H^2 * |A| * (1 - exp(-lambda)).

    # The instruction "output each number in the final equation" is interpreted
    # as clearly presenting the components of the final symbolic formula.
    # The numbers '2' (from the exponent) and '1' are inherent parts of the expression.
    final_expression = f"{H}^2 * {A_size} * (1 - exp(-{lam}))"

    print("The tightest upper bound for J(pi^*) - J(pi_hat) is given by the expression:")
    print(final_expression)

# Execute the function to print the solution.
solve_imitation_learning_bound()