import math

def sum_prime_factors(n):
    """
    Calculates the sum of the prime factors of a number n, with multiplicity.
    For example, sum_prime_factors(12) = 2 + 2 + 3 = 7.
    """
    s = 0
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            s += d
            temp_n //= d
        d += 1
    if temp_n > 1:
        s += temp_n
    return s

def solve():
    """
    Finds the number N based on the logic derived from the problem description.
    """
    # Step 1: Create a map from a sum-of-prime-factors value (s) to the smallest
    # integer (n) that has that sum. We iterate n upwards to ensure we find the smallest.
    # A limit of 100 for n should be sufficient to find the small numbers required.
    spf_to_min_n = {}
    limit = 100 
    for n in range(2, limit):
        s = sum_prime_factors(n)
        if s not in spf_to_min_n:
            spf_to_min_n[s] = n

    # Step 2: Find the pair (a, b) that satisfies the conditions and minimizes a + b.
    min_sum = float('inf')
    result_a = -1
    result_b = -1

    # Iterate through all possible sums s_a that a can have.
    # We only need to check s_a up to 10, as s_b = 20 - s_a would be symmetric.
    for s_a in range(2, 11): 
        s_b = 20 - s_a
        
        # Check if we have found numbers for both s_a and s_b
        if s_a in spf_to_min_n and s_b in spf_to_min_n:
            a = spf_to_min_n[s_a]
            b = spf_to_min_n[s_b]
            
            # Condition ¬R: a and b must be different.
            if a == b:
                continue

            current_sum = a + b
            if current_sum < min_sum:
                min_sum = current_sum
                result_a = a
                result_b = b

    # Step 3: Print the result.
    # The number N is the smallest sum found.
    # The problem asks for the numbers in the final equation.
    if result_a != -1:
        # Sort to present the equation in a standard way (e.g., 3 + 17, not 17 + 3)
        final_a = min(result_a, result_b)
        final_b = max(result_a, result_b)
        n_val = final_a + final_b
        print(f"The number N is the sum in the equation: {final_a} + {final_b} = {n_val}")
    else:
        print("No solution found within the checked limit.")

solve()