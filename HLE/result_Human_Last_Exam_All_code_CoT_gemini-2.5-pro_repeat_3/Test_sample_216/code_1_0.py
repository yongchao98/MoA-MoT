def solve():
    """
    This function formulates and prints the tightest upper bound for the performance difference
    J(pi^*) - J(hat{pi}) based on the provided information.

    The steps are as follows:
    1.  The problem setting involves an imitation learning algorithm where the risk is measured on the learner's state visitation distribution. This is characteristic of algorithms like DAgger, which have performance bounds that scale linearly with the horizon H.
    2.  The standard performance difference bound for such algorithms is J(pi^*) - J(hat{pi}) <= H * T(hat{pi}, pi^*), where T is the population TV risk and we assume rewards are normalized to have a maximum value of 1.
    3.  The problem provides a specific bound for this risk: T(hat{pi}, pi^*) <= |A| * (1 - exp(-lambda)).
    4.  By substituting the given risk bound into the performance difference bound, we obtain the final expression.
    """

    # Define the symbols for the final equation
    H = "H"
    A_size = "|A|"
    lmbda = "lambda"
    one = "1"
    
    # The final bound is the product of the horizon H and the given upper bound on the TV risk.
    # The equation is H * |A| * (1 - exp(-lambda)).
    # We print each part of the formula as requested.
    print(f"The tightest upper bound for J(pi^*) - J(hat{pi}) is:")
    print(f"{H} * {A_size} * ({one} - exp(-{lmbda}))")

solve()