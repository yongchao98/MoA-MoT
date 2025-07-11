import math
import sys

# Set a higher recursion limit for large m
sys.setrecursionlimit(2000)

memo_c = {}

def combinations(n, k, cache):
    """Helper function for combinations to handle negative k"""
    if k < 0 or k > n:
        return 0
    if (n, k) in cache:
        return cache[(n, k)]
    res = math.comb(n, k)
    cache[(n, k)] = res
    return res

def calculate_c(m, A):
    """
    Recursively calculates the auxiliary term c_m.
    A = 2**(n-1) - 1
    """
    if m in memo_c:
        return memo_c[m]
    
    if m == 0:
        return 1
    if m == 1:
        return -1

    c_prev = calculate_c(m - 1, A)
    
    if m % 2 == 1: # m is odd
        result = -c_prev
    else: # m is even
        k = m // 2
        # This division is always exact.
        result = (c_prev * (A - k + 1)) // k
        
    memo_c[m] = result
    return result

def solve_sets_problem(n, m):
    """
    Calculates the number of sets T satisfying the conditions.
    """
    # According to the problem statement, n and m are positive integers.
    if n <= 0 or m <= 0:
        print("n and m must be positive integers.")
        return

    # Total number of non-empty subsets of S
    num_subsets = 2**n - 1

    if m > num_subsets:
        # Cannot choose m distinct subsets if m is too large
        print(0)
        return

    # Calculate the main terms of the formula
    comb_term = combinations(num_subsets, m, {})

    # The parameter A for the c_m calculation
    A = 2**(n-1) - 1
    
    # Clear memoization cache for c_m calculation
    global memo_c
    memo_c = {}
    c_m = calculate_c(m, A)
    
    c_term = num_subsets * c_m

    # The final formula ensures exact division
    numerator = comb_term + c_term
    denominator = 2**n
    
    result = numerator // denominator
    
    # Print the full equation for clarity
    print(f"f(n={n}, m={m}) = ( C({num_subsets}, {m}) + ({num_subsets}) * ({c_m}) ) / {denominator}")
    print(f"f({n}, {m}) = ( {comb_term} + {c_term} ) / {denominator}")
    print(f"f({n}, {m}) = {numerator} / {denominator}")
    print(f"Result: {result}")


if __name__ == '__main__':
    # Example from the problem description, change these values as needed
    # Positive integers n and m
    n_val = 4
    m_val = 3
    
    solve_sets_problem(n_val, m_val)

    # Another example: n=3, m=4
    # n_val = 3
    # m_val = 4
    # solve_sets_problem(n_val, m_val)