import math

def sum_prime_factors(n):
    """
    Calculates the sum of the prime factors of a number n, including multiplicity.
    For example, sum_prime_factors(12) = 2 + 2 + 3 = 7.
    """
    if n <= 1:
        return 0
    
    total_sum = 0
    temp_n = n
    
    # Handle factor 2
    while temp_n % 2 == 0:
        total_sum += 2
        temp_n //= 2
        
    # Handle odd factors
    d = 3
    while d * d <= temp_n:
        while temp_n % d == 0:
            total_sum += d
            temp_n //= d
        d += 2
        
    # If temp_n is still greater than 1, it must be a prime factor
    if temp_n > 1:
        total_sum += temp_n
        
    return int(total_sum)

def find_smallest_sum():
    """
    Finds the smallest integer N which is the sum of two different integers a and b,
    where the sum of the prime factors of a and b is 20.
    """
    min_sum = float('inf')
    result_a = -1
    result_b = -1

    # We search for a and b in a reasonable range. 
    # A small N implies small a and b, so 100 is a safe upper bound.
    search_limit = 100

    # Cache results for sum_prime_factors to avoid re-computation
    spf_cache = {i: sum_prime_factors(i) for i in range(1, search_limit)}

    for a in range(1, search_limit):
        # We need b to be different from a. Looping from a + 1 ensures this
        # and avoids checking pairs twice (e.g., (3,17) and (17,3)).
        for b in range(a + 1, search_limit):
            spf_a = spf_cache[a]
            spf_b = spf_cache[b]

            # Condition from ¬Q: sum of prime factors of a and b is 20
            if spf_a + spf_b == 20:
                current_sum = a + b
                # Condition from ¬P: N is the *smallest* such sum
                if current_sum < min_sum:
                    min_sum = current_sum
                    result_a = a
                    result_b = b
    
    # Print the final equation as requested
    print(f"{result_a} + {result_b} = {min_sum}")

find_smallest_sum()
<<<20>>>