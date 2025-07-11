import math

def get_sum_of_distinct_prime_factors(n):
    """
    Calculates the sum of distinct prime factors of a given integer n.
    For example:
    - get_sum_of_distinct_prime_factors(12) = 2 + 3 = 5
    - get_sum_of_distinct_prime_factors(1) = 0
    """
    if n <= 1:
        return 0
    
    factors = set()
    d = 2
    temp_n = n
    
    # Check for factors up to the square root of n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    
    # If there's a remainder, it's a prime factor itself
    if temp_n > 1:
        factors.add(temp_n)
        
    return sum(factors)

def find_smallest_n_and_print_equation():
    """
    Finds the smallest integer N that is the sum of two different integers a and b,
    where the sum of the distinct prime factors of a and b is 20.
    """
    # We search for the sum N, starting from the smallest possible sum of two different positive integers (1+2=3)
    n_sum = 3
    while True:
        # For a given sum N, check all pairs (a, b) where a < b and a + b = N.
        # We only need to iterate 'a' up to half of N to cover all pairs.
        for a in range(1, (n_sum // 2) + 1):
            b = n_sum - a

            # Ensure a and b are different. This is guaranteed if a != n_sum/2.
            if a == b:
                continue
            
            sopf_a = get_sum_of_distinct_prime_factors(a)
            sopf_b = get_sum_of_distinct_prime_factors(b)
            
            # Check if the sum of prime factors condition is met
            if sopf_a + sopf_b == 20:
                # Since we are iterating n_sum upwards, the first match is the smallest N.
                print(f"Found the smallest number N = {n_sum}.")
                print(f"It is the sum of a = {a} and b = {b}.")
                print(f"The sum of the distinct prime factors of {a} is {sopf_a}.")
                print(f"The sum of the distinct prime factors of {b} is {sopf_b}.")
                print(f"The total sum of prime factors is {sopf_a} + {sopf_b} = 20.")
                print("\nThe final equation is:")
                print(f"{a} + {b} = {n_sum}")
                return

# Execute the function to find the number and print the result.
find_smallest_n_and_print_equation()
<<<20>>>