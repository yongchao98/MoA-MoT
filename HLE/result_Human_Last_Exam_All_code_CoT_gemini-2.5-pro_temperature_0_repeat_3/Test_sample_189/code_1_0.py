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

def find_prime_chain():
    """
    Finds the largest prime p < 1,000,000 such that p = 4u+1 and u = 4v+1,
    where u and v are also prime.
    """
    # The computer's memory limit for an integer is 999,999.
    limit = 1000000
    
    # From p = 16v+5, the maximum possible value for v is (limit-5)/16.
    max_v = (limit - 5) // 16
    
    largest_p = 0
    corresponding_u = 0
    corresponding_v = 0

    # We iterate through all possible values for v to find all valid chains.
    for v in range(2, max_v + 1):
        # Condition 1: v must be prime.
        if is_prime(v):
            u = 4 * v + 1
            # Condition 2: u must be prime.
            if is_prime(u):
                p = 4 * u + 1
                # Condition 3: p must be prime and within the computer's limit.
                if p < limit and is_prime(p):
                    # Since we are iterating upwards, the last found values will be the largest.
                    largest_p = p
                    corresponding_u = u
                    corresponding_v = v
    
    # Print the final result in the format p:u:v
    # The problem asks to output each number in the final equation.
    print(f"{largest_p}:{corresponding_u}:{corresponding_v}")

find_prime_chain()