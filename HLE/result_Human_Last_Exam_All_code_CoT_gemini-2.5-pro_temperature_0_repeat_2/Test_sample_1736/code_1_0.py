import math

# Memoization cache to store results and avoid re-calculating.
memo = {}

def sum_prime_factors(n):
    """
    Calculates the sum of the prime factors of a number n, with multiplicity.
    For example, sum_prime_factors(12) = 2 + 2 + 3 = 7.
    Uses memoization to speed up repeated calculations for the same number.
    """
    if n in memo:
        return memo[n]
    
    if n <= 1:
        return 0

    original_n = n
    total_sum = 0
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            total_sum += d
            temp_n //= d
        d += 1
    if temp_n > 1:
        total_sum += temp_n
    
    memo[original_n] = total_sum
    return total_sum

def find_smallest_sum():
    """
    This function finds the number N that satisfies the given logical conditions.
    It searches for the smallest integer N which is a sum of two different integers
    a and b, where the sum of the prime factors of a and b is 20.
    """
    target_prime_factor_sum = 20
    
    # We iterate through possible sums N, starting from 3.
    # The upper limit of 100 is sufficient as the answer is expected to be small.
    for n_sum in range(3, 100): 
        # For each sum N, iterate through possible values of 'a'.
        # 'a' goes up to n_sum / 2, and b is the remainder.
        # This ensures a < b and we don't check pairs twice (e.g., (3,17) and (17,3)).
        for a in range(1, n_sum // 2 + 1):
            b = n_sum - a
            
            # The loop structure a < n_sum/2 ensures a < b, so they are different.
            if a == b:
                continue

            # Check if the sum of prime factors equals our target.
            spf_a = sum_prime_factors(a)
            spf_b = sum_prime_factors(b)
            
            if spf_a + spf_b == target_prime_factor_sum:
                # Since we are iterating n_sum upwards, the first one we find is the smallest.
                print(f"Found the smallest number N = a + b satisfying the conditions.")
                print(f"The two different integers are a = {a} and b = {b}.")
                print(f"The sum of prime factors for a={a} is {spf_a}.")
                print(f"The sum of prime factors for b={b} is {spf_b}.")
                print(f"The total sum of prime factors is {spf_a} + {spf_b} = {target_prime_factor_sum}.")
                print("The final equation is:")
                print(f"{a} + {b} = {n_sum}")
                return n_sum
    return None

# Execute the function to find the number.
result = find_smallest_sum()
