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

def find_special_prime_chain():
    """
    Finds the largest prime p < 1,000,000 such that:
    p = 4u + 1
    u = 4v + 1
    and p, u, v are all prime.
    """
    # The memory size of 999999 suggests a search limit for p.
    limit_p = 1000000
    
    # From p = 16v + 5, we can derive the limit for v.
    limit_v = (limit_p - 5) // 16

    # We search downwards from the highest possible v to find the largest p first.
    for v in range(limit_v, 1, -1):
        # Check 1: Is v prime?
        if is_prime(v):
            u = 4 * v + 1
            # Check 2: Is u = 4v + 1 prime?
            if is_prime(u):
                p = 4 * u + 1
                # Check 3: Is p = 4u + 1 prime and within the limit?
                if p < limit_p and is_prime(p):
                    # Since we are searching downwards, the first result found is the largest.
                    # The problem requires printing each number in the final equation.
                    print(f"{p}:{u}:{v}")
                    return

    print("No such prime chain found within the specified limit.")

if __name__ == '__main__':
    find_special_prime_chain()