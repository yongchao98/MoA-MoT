def solve_complexity():
    """
    Calculates and prints the complexity tuples for the two regimes.
    """
    # For Regime 1: N = 2^sqrt(L), the complexity is Theta(N log N)
    # (N log N)^2 = N^2 (log N)^2
    a1 = 2
    b1 = 2
    c1 = 0

    # For Regime 2: N = 2^((log L)^2), the complexity is also Theta(N log N)
    # (N log N)^2 = N^2 (log N)^2
    a2 = 2
    b2 = 2
    c2 = 0

    # The format is (a,b,c),(a,b,c)
    result = f"({a1},{b1},{c1}),({a2},{b2},{c2})"
    print(result)

solve_complexity()