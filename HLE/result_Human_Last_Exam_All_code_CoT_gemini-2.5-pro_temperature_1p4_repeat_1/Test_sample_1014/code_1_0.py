import math

def find_largest_prime_factor(n):
    """
    Finds the largest prime factor of a given integer n and prints the
    prime factorization equation.
    """
    original_n = n
    factors = []

    # Handle the factor 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors
    # We iterate from 3 up to the square root of the remaining n
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # This condition is to handle the case when n is a prime number
    # greater than 2, which would be the last remaining factor.
    if n > 1:
        factors.append(n)

    # Building the equation string from the factors
    # e.g., 8139881 = 131 * 62137
    equation = f"{original_n} = {' * '.join(map(str, factors))}"
    print("The prime factorization is:")
    print(equation)

    if factors:
        largest_prime = max(factors)
        print(f"The largest prime divisor of {original_n} is: {largest_prime}")
    else:
        # This case handles if the original number is 1 or 0
        print(f"{original_n} has no prime divisors.")

# The number we want to factor
number_to_factor = 8139881
find_largest_prime_factor(number_to_factor)