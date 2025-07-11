def solve():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n.
    """
    n = 35

    # The formula for the minimum sum of vertices is n * (n + 1) / 2 - 1.
    # This is derived from the optimal "star decomposition" of K_n.
    numerator = n * (n + 1)
    denominator = 2
    
    # We perform the calculation step-by-step to show the process.
    result_before_subtraction = numerator // denominator
    final_result = result_before_subtraction - 1
    
    print(f"The problem is to find the minimum number of vertices in a family of bipartite graphs covering the edges of K_n for n = {n}.")
    print("This is a known result in graph theory, and the minimum value is achieved by a star decomposition.")
    print(f"The formula for this minimum is: n * (n + 1) / 2 - 1")
    print("\nCalculation steps:")
    print(f"({n} * ({n} + 1)) / {denominator} - 1")
    print(f"= ({n} * {n + 1}) / {denominator} - 1")
    print(f"= {numerator} / {denominator} - 1")
    print(f"= {result_before_subtraction} - 1")
    print(f"= {final_result}")

solve()