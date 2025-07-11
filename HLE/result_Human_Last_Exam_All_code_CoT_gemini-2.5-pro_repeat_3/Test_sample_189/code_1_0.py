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
    # Loop from 5 up to the square root of n, with a step of 6
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_special_primes():
    """
    Finds the largest prime p of the form p = 4u+1, where u is a prime
    of the form u = 4v+1, and p is within the computer's memory limit.
    """
    max_p_found = -1
    corresponding_u = -1
    corresponding_v = -1

    # The largest number p cannot exceed 999999.
    # p = 4u + 1 = 4(4v + 1) + 1 = 16v + 5
    # 16v + 5 <= 999999  =>  16v <= 999994  =>  v <= 62499.625
    limit_v = 62499

    # Iterate v from the largest possible value downwards to find the max p first.
    for v in range(limit_v, 0, -1):
        # Step 1: Check if v is a prime.
        if is_prime(v):
            u = 4 * v + 1
            # u = 4v+1 is already in the form 4k+1.
            # Step 2: Check if u is a prime.
            if is_prime(u):
                p = 4 * u + 1
                # p = 4u+1 is already in the form 4k+1.
                # The loop for v ensures p will be <= 999999.
                # Step 3: Check if p is a prime.
                if is_prime(p):
                    # Since we are iterating downwards, the first one we find is the largest.
                    max_p_found = p
                    corresponding_u = u
                    corresponding_v = v
                    # We can stop now because we found the largest possible p.
                    break
    
    # Print the final result in the specified format p:u:v
    if max_p_found != -1:
        print(f"The largest prime p is {max_p_found}, derived from u = {corresponding_u} and v = {corresponding_v}.")
        print("Final result in p:u:v format:")
        print(f"{max_p_found}:{corresponding_u}:{corresponding_v}")
    else:
        print("No such prime triplet was found within the given constraints.")

# Execute the function to find and print the result.
find_special_primes()