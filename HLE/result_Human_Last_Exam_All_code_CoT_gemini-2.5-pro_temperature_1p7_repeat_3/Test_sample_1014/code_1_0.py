import math

def solve():
    """
    Finds the largest prime divisor of 8139881.
    """
    n = 8139881
    original_n = n
    factors = []

    # Handle the factor 2, if any
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors from 3 upwards
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is still greater than 1 at this point,
    # the remaining n must be a prime factor itself.
    if n > 1:
        factors.append(n)

    # To satisfy the request to "output each number in the final equation",
    # we print the prime factorization.
    # If the number is prime, the equation will be "n = n".
    equation_str = f"{original_n} = {' * '.join(map(str, factors))}"
    print(equation_str)

    # Get the largest factor from the list
    if factors:
        largest_prime = max(factors)
        print(f"The largest prime divisor of {original_n} is {largest_prime}.")
    else:
        # This case applies to numbers <= 1
        print(f"{original_n} has no prime divisors.")

solve()
<<<8139881>>>