import math

def solve_imitation_learning_bound():
    """
    This function derives and prints the tightest upper bound for the performance
    difference between an expert and a learned policy in imitation learning.

    The derivation follows these steps:
    1. The performance difference J(pi^*) - J(pi_hat) is bounded by the
       sum of accumulated one-step deviation errors.
    2. A one-step error at time t is bounded by (H-t) * TV_risk(t), assuming
       rewards are in [0, 1].
    3. The total error is sum_{t=0 to H-1} of (H-t) * E[TV_risk(t)].
    4. To get the tightest upper bound given a constraint on the average TV risk,
       we consider the worst-case where the error is concentrated at t=0.
       This leads to a dependency of H^2 * (Average TV Risk).
    5. The problem states that the population TV risk is bounded by
       |A|*(1 - e^-lambda).
    """

    # Symbolic representation of the variables
    H = "H"
    abs_A = "|A|"
    lam = "lambda"

    # Constructing the final expression for the tightest upper bound.
    # The exponent is 2.
    exponent = "2"
    # The constant for the first term in the parentheses is 1.
    constant_one = "1"
    # The negative sign for lambda
    neg_one = "-1" # conceptually, but just displayed as "-"

    # The formula is H^2 * population_tv_risk_bound
    # H^2 comes from the error accumulation in imitation learning.
    # The bound for population TV risk is given.
    bound_formula = f"{H}**{exponent} * {abs_A} * ({constant_one} - exp(-{lam}))"

    print("The tightest upper bound of J(pi^*) - J(pi_hat) is given by the expression:")
    print(bound_formula)
    print("\nWhere:")
    print(f"  {H}: the horizon of the episode")
    print(f"  {abs_A}: the size of the discrete action space")
    print(f"  {lam}: a hyperparameter of the algorithm")
    print(f"  exp: the exponential function")

solve_imitation_learning_bound()