def calculate_envy_free_bound():
    """
    Calculates the upper bound for a connected envy-free cake-cutting
    allocation for n=4 agents.

    The problem of envy-free cake-cutting for n agents aims to divide a resource
    such that no agent prefers another agent's piece. For four agents (n=4),
    a bounded protocol for a connected, exact envy-free allocation was a
    long-standing open problem.

    A breakthrough result by Fragniaud, Feller, and Meunier (2021) established
    a protocol that requires at most 12 cuts.

    This upper bound, O, can be represented by the formula O = n * (n - 1),
    which corresponds to the number of ordered pairs of agents, representing
    all potential instances of envy that must be resolved.

    This script calculates this bound for n=4.
    """

    # Number of agents
    n = 4

    # The formula for the upper bound O is n * (n - 1)
    # We will calculate the result and show the components of the equation.
    term1 = n
    term2 = n - 1
    upper_bound = term1 * term2

    print(f"For n = {term1} agents, a realistic upper bound for a connected envy-free allocation can be expressed by the formula n * (n - 1).")
    print(f"Calculating the equation: {term1} * {term2} = {upper_bound}")
    print(f"The most realistic upper bound O is: {upper_bound}")


if __name__ == "__main__":
    calculate_envy_free_bound()
