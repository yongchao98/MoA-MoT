def solve_cake_cutting_bound():
    """
    Calculates an upper bound for the envy-free cake-cutting problem
    for n=4 based on an extrapolated formula.
    """
    # The number of agents in the problem.
    n = 4

    # The formula is based on the pattern n^(n-1) observed for n=2 and n=3.
    exponent = n - 1
    result = n ** exponent

    # Print the equation step-by-step as requested.
    print(f"Based on the pattern n^(n-1), the upper bound for n = {n} agents can be calculated.")
    print(f"The equation is: {n}^({n} - 1)")
    print(f"Which simplifies to: {n}^{exponent}")
    print(f"The final result is: {result}")

solve_cake_cutting_bound()