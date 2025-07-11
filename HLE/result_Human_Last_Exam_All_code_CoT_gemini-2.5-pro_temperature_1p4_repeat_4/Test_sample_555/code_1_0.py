import sys
from fractions import Fraction

# It's possible for the recursion to go deep, so increase the limit.
sys.setrecursionlimit(2000)

# Memoization cache to store results for E(n)
# The value will be a tuple: (expected_rolls, best_k, remainder_r)
memo = {}

def get_expected_rolls(n):
    """
    Calculates the minimum expected number of rolls E(n) to select one of n items
    uniformly using a 7-sided die.
    """
    if n in memo:
        return memo[n][0]

    # Base case: If there's 1 or 0 choice, no rolls are needed.
    if n <= 1:
        return Fraction(0)

    min_expected_value = float('inf')
    best_k = -1
    best_r = -1
    
    # We need to find the best k that minimizes the expected value.
    # Let's check k from 1 up to a reasonable limit (e.g., 10).
    # The term 'k' will dominate, so we don't need a large k.
    for k in range(1, 10):
        d_k = 7**k
        if d_k < n:
            continue
        
        r = d_k % n
        
        # E(n) = k + (r / d^k) * E(r)
        current_expected_value = Fraction(k) + Fraction(r, d_k) * get_expected_rolls(r)
        
        if current_expected_value < min_expected_value:
            min_expected_value = current_expected_value
            best_k = k
            best_r = r

    # Store the result in the cache
    memo[n] = (min_expected_value, best_k, best_r)
    return min_expected_value

def solve():
    """
    Solves the problem and prints the steps and final answer.
    """
    # Calculate E(10)
    final_result = get_expected_rolls(10)
    
    print("The optimal strategy is found by solving a recursive equation.")
    print("Let E(n) be the minimal expected number of rolls to choose from n items.")
    print("The derivation for E(10) is as follows:\n")

    path = []
    curr_n = 10
    while curr_n > 1:
        if curr_n in memo:
            res, k, r = memo[curr_n]
            d_k = 7**k
            path.append(f"E({curr_n}) = {k} + ({r} / {d_k}) * E({r})")
            if r == 0:
                path.append(f"E(0) = 0, so the expression simplifies.")
            elif r == 1:
                 path.append(f"E(1) = 0, so the expression simplifies.")
            curr_n = r
        else:
            break
            
    # Reverse the path to show it from the bottom up for clarity, although top-down is also fine.
    # For this problem, let's keep it top-down.
    for line in path:
        print(line)

    print(f"\nBy solving this system, we get the final value for E(10).")
    print(f"The minimal expected value of rolls is E(10) = {final_result.numerator}/{final_result.denominator}.")


solve()