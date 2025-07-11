import math

def is_prime(n):
    """
    Checks if a number n is prime using an efficient trial division method.
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
    Finds the largest prime p < 1,000,000 of the form p=4u+1, where u is
    a prime of the form u=4v+1.
    """
    # Based on p < 1,000,000, the maximum possible value for v is derived:
    # p = 4u + 1 < 1000000 => 4u < 999999 => u < 249999.75
    # u = 4v + 1 < 250000 => 4v < 249999 => v < 62499.75
    # So we start v from 62499 and search downwards.
    v_max = 62499

    for v in range(v_max, 0, -1):
        u = 4 * v + 1
        
        # Check if u is prime
        if is_prime(u):
            p = 4 * u + 1
            
            # The value of p is guaranteed to be < 1,000,000 by the v_max calculation
            # but we check just in case and to ensure it's prime.
            if is_prime(p):
                # Since we are iterating downwards, the first result we find is the largest.
                print(f"{p}:{u}:{v}")
                return

if __name__ == "__main__":
    find_largest_special_prime()
