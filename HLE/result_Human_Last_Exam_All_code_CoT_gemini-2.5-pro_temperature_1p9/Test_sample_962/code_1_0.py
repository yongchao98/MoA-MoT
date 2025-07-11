import math

def print_prime_factorization(n):
    """
    Calculates the prime factorization of a given number n and prints it
    in the format n = f1 * f2 * ...
    """
    if not isinstance(n, int) or n <= 1:
        print("Please provide an integer greater than 1.")
        return

    original_n = n
    factors = []

    # Handle factors of 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors
    # We check from 3 up to the square root of the remaining number n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n //= i

    # This condition is to handle the case when n is a prime number
    # greater than 2 at the end of the loop
    if n > 2:
        factors.append(n)

    # Convert factors to strings for printing the equation
    str_factors = [str(f) for f in factors]
    
    # Print the final equation
    print(f"{original_n} = {' * '.join(str_factors)}")

# Example Usage with the number 360
number_to_factorize = 360
print_prime_factorization(number_to_factorize)