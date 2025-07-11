def solve():
    """
    This function calculates and prints the tightest upper bound for the performance
    difference J(pi^*) - J(hat_pi) based on the provided information.

    The derivation is as follows:
    1. A standard result in imitation learning bounds the performance difference due to
       the compounding error effect. This bound relates the difference in expected
       returns to the total variation (TV) risk between the expert and learned policies.
       Assuming rewards are normalized to be within [0, 1] (i.e., R_max = 1), a common
       form of this bound is:
       J(pi^*) - J(hat_pi) <= 2 * H**2 * T(hat_pi, pi^*)
       where T(hat_pi, pi^*) is the population TV risk, defined as the expected TV
       distance over the state distribution induced by the expert policy pi^*.

    2. The problem provides an upper bound on this TV risk, which depends on a
       hyperparameter lambda and the action space size |A|:
       T(hat_pi, pi^*) <= |A| * (1 - exp(-lambda))

    3. By substituting the second inequality into the first, we get the tightest
       upper bound on the performance difference given the available information:
       J(pi^*) - J(hat_pi) <= 2 * H**2 * |A| * (1 - exp(-lambda))

    The code below will print this final expression.
    """
    # Define the symbols used in the expression for clarity.
    # H: Horizon of the episode
    # |A|: Size of the discrete action space
    # lambda: A hyperparameter of the algorithm
    # exp: The exponential function
    
    # The final expression for the upper bound.
    # Note: '**' denotes exponentiation.
    upper_bound_expression = "2 * H**2 * |A| * (1 - exp(-lambda))"
    
    print("The tightest upper bound of J(pi^*) - J(hat_pi) is given by the expression:")
    print(upper_bound_expression)

solve()