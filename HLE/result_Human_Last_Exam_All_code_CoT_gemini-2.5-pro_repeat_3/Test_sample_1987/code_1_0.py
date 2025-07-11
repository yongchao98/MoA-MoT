import math
from fractions import Fraction

def calculate_f(args):
    """
    Calculates the value of f(a_1, ..., a_n) using the hook-content formula.
    The arguments are expected to be sorted: a_1 <= a_2 <= ... <= a_n.
    """
    # Check for the special case where all arguments are equal.
    # In this case, the Young diagram is a rectangle, and the formula simplifies to 1.
    if len(set(args)) == 1:
        return 1
        
    n = len(args)
    # The partition lambda is the reverse of the sorted arguments.
    # The input 'args' is already sorted.
    lambda_part = list(reversed(args))
    
    # Calculate the transpose of the partition (lambda'), which gives the column lengths.
    if not lambda_part:
        return 1 # f() = 1 based on f(0...0)=1 for n=0
    max_col = lambda_part[0]
    lambda_prime = [0] * max_col
    for part_size in lambda_part:
        for i in range(part_size):
            lambda_prime[i] += 1
            
    total_value = Fraction(1)
    
    # Iterate through each cell (r, c) of the Young diagram for lambda.
    # r and c are 1-indexed.
    for r in range(1, len(lambda_part) + 1):
        for c in range(1, lambda_part[r-1] + 1):
            # content(u) = c - r
            content = c - r
            
            # hook_length(u) = (lambda_r - c) + (lambda'_c - r) + 1
            hook_length = (lambda_part[r-1] - c) + (lambda_prime[c-1] - r) + 1
            
            # Term in the product is (n + content) / hook_length
            term = Fraction(n + content, hook_length)
            total_value *= term
            
    # The result must be an integer.
    if total_value.denominator != 1:
        # This should not happen for valid inputs.
        raise ValueError("Calculation resulted in a non-integer value.")
        
    return total_value.numerator

def solve():
    """
    Solves the user's request by calculating the three specified values.
    """
    # 1. Calculate f(2, 4, 5)
    ans1 = calculate_f((2, 4, 5))
    
    # 2. Calculate f(9000, 9000, 9000)
    # This falls into the special case f(A, A, A) = 1
    ans2 = calculate_f((9000, 9000, 9000))
    
    # 3. Calculate f(p, p, p, p) mod p for p = 10^9 + 7
    # This is also a special case f(A, A, A, A) = 1.
    # So, 1 mod p = 1.
    p = 10**9 + 7
    val_f_p = calculate_f((p, p, p, p))
    ans3 = val_f_p % p
    
    print(f"{ans1},{ans2},{ans3}")

solve()