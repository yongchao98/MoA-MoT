import math

def find_largest_prime_divisor():
    """
    Finds all prime factors of 8139881, prints the factorization,
    and identifies the largest prime factor.
    """
    n = 8139881
    original_n = n
    factors = []

    # Handle the factor 2, if any
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors starting from 3
    # We only need to check up to the square root of the current n
    i = 3
    while i * i <= n:
        if n % i == 0:
            factors.append(i)
            n //= i
        else:
            i += 2

    # If n is still greater than 1 at the end, it is a prime factor itself.
    if n > 1:
        factors.append(n)

    # Determine the result and print the output
    if not factors:
        # This case would only occur for n=1
        largest_prime = 1
        print(f"The number {original_n} has no prime factors.")
    elif len(factors) == 1:
        # This is the case for a prime number
        largest_prime = original_n
        print(f"{original_n} is a prime number.")
    else:
        # Print the full prime factorization equation
        equation = f"{original_n} = " + " * ".join(map(str, factors))
        print(equation)
        largest_prime = max(factors)

    print(f"\nThe largest prime divisor of {original_n} is {largest_prime}.")


# Run the function to solve the problem
find_largest_prime_divisor()

<<<8139881>>>