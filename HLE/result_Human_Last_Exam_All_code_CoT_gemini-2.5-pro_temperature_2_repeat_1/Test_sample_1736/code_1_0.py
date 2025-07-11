import math

def get_sum_of_prime_factors(n):
    """
    Calculates the sum of the prime factors of a number n, with multiplicity.
    For example, for n=12 (2*2*3), the sum is 2+2+3=7.
    Returns 0 for n <= 1.
    """
    if n <= 1:
        return 0
    
    total_sum = 0
    temp_n = n
    
    # Handle factor of 2
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
        
    # This case happens if temp_n is a prime number greater than 2
    if temp_n > 1:
        total_sum += temp_n
        
    return total_sum

def find_smallest_N():
    """
    Finds the smallest integer N that is the sum of two different integers a and b,
    where the sum of the prime factors of a and b is 20.
    """
    # Iterate through possible values of N, starting from the smallest sums.
    # An upper bound of 100 is more than sufficient.
    for N in range(2, 100):
        # To find a and b, iterate a from 1 up to N/2.
        # This ensures a < b and avoids duplicate pairs.
        for a in range(1, math.ceil(N / 2)):
            b = N - a

            # This condition is already guaranteed by the loop range `range(1, ceil(N/2))`.
            if a == b:
                continue

            # Check if the sum of prime factors is 20.
            if get_sum_of_prime_factors(a) + get_sum_of_prime_factors(b) == 20:
                # Since we iterate N from small to large, the first one found is the smallest.
                # Output the numbers in the final equation as requested.
                print(f"{N} = {a} + {b}")
                print(f"<<<{N}>>>")
                return

find_smallest_N()