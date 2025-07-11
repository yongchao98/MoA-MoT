import sys

# A cache to store the results of sum_prime_factors to speed up the process
memo = {}

def sum_prime_factors(n):
    """
    Calculates the sum of the prime factors of a number n, with multiplicity.
    For example, sum_prime_factors(12) = 2 + 2 + 3 = 7.
    Uses a cache (memo) to store results.
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
    if temp_n > 1:  # The remainder is a prime factor
        total_sum += temp_n

    memo[original_n] = total_sum
    return total_sum

def find_smallest_n_with_factor_sum(target_sum):
    """
    Finds the smallest integer N that is the sum of two different integers a and b,
    such that the sum of the prime factors of a and b equals target_sum.
    """
    n = 3  # Smallest possible sum of two different positive integers (1+2)
    while True:
        # Iterate through a from 1 up to n/2.
        # This ensures a < b, so a and b are different.
        for a in range(1, (n // 2) + 1):
            b = n - a
            if a == b:
                continue

            s_a = sum_prime_factors(a)
            s_b = sum_prime_factors(b)

            if s_a + s_b == target_sum:
                print(f"Found the smallest N satisfying the conditions.")
                print(f"a = {a}, b = {b}")
                print(f"Sum of prime factors of {a} is: {s_a}")
                print(f"Sum of prime factors of {b} is: {s_b}")
                print(f"Total sum of prime factors: {s_a} + {s_b} = {s_a + s_b}")
                print("\nFinal equation that gives the smallest N:")
                print(f"{a} + {b} = {n}")
                return n
        n += 1

# The condition from ¬Q is that the sum of prime factors is 20.
TARGET_FACTOR_SUM = 20
result_n = find_smallest_n_with_factor_sum(TARGET_FACTOR_SUM)
print(f"\nThe number N that satisfies the negation of P v R v Q is {result_n}.")
