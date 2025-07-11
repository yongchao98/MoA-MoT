def solve_ratio_limit():
    """
    Calculates the expected number of remaining elements and the ratio to n.

    Let E_n be the expected number of remaining elements for a row of size n.
    The recurrence relation is derived as:
    E_n = ((n - 2) * E_{n-1} + 2 * E_{n-2}) / (n - 1) for n >= 2
    with base cases E_0 = 0 and E_1 = 1.

    This function calculates E_n/n for several large n to observe the limit.
    """
    
    # We use a list to store the values of E_n to avoid recomputation (memoization).
    # Using an iterative approach to avoid Python's recursion depth limit.
    e = [0.0] * 20001
    e[0] = 0
    e[1] = 1

    # Print header
    print(f"{'n':>8} | {'E_n':>15} | {'E_n / n':>15}")
    print("-" * 43)

    # Calculate E_n iteratively
    for n in range(2, 20001):
        e[n] = ((n - 2) * e[n - 1] + 2 * e[n - 2]) / (n - 1)
        
        # Print the ratio for specific values of n to show the trend
        if n in [10, 100, 1000, 5000, 10000, 20000]:
            ratio = e[n] / n
            print(f"{n:>8} | {e[n]:>15.4f} | {ratio:>15.8f}")

    final_ratio = e[20000] / 20000
    print("\nAs n approaches infinity, the ratio E_n / n appears to approach 0.")
    # The final equation asked for is the limit value
    print(f"The limit of the ratio is: {0}")


solve_ratio_limit()