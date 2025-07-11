import math

def get_sum_of_prime_factors(n):
    """
    Calculates the sum of the prime factors of a given number n, including multiplicity.
    For example, for 12 = 2*2*3, the sum is 2+2+3=7.
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
    for i in range(3, int(math.sqrt(temp_n)) + 1, 2):
        while temp_n % i == 0:
            total_sum += i
            temp_n //= i
            
    # If temp_n is still a prime number greater than 2
    if temp_n > 2:
        total_sum += temp_n
        
    return total_sum

def find_smallest_number():
    """
    Finds the smallest integer N that is the sum of two distinct integers a and b,
    such that the sum of their prime factors is 20.
    """
    n = 1
    while True:
        # We check for a and b such that a + b = n and a < b
        # We only need to check a up to n/2
        for a in range(1, n // 2 + 1):
            b = n - a
            
            # Condition: a and b must be different.
            # By iterating a up to n/2 and setting b = n-a, a will always be <= b.
            # We skip the case where a == b.
            if a == b:
                continue

            # Condition: sum of prime factors of a and b is 20
            if get_sum_of_prime_factors(a) + get_sum_of_prime_factors(b) == 20:
                # We found the smallest N, along with the corresponding a and b.
                # Print the final equation as requested.
                print(f"{n} = {a} + {b}")
                return
        n += 1

find_smallest_number()