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

def find_largest_prime_chain():
    """
    Finds the largest prime p < 1,000,000 of the form p = 4u + 1,
    where u is a prime of the form u = 4v + 1, and v is also a prime.
    """
    # From p = 16v + 5 and p < 1,000,000, we find v_max.
    v_max = (1000000 - 5) // 16

    # We search downwards from v_max to find the largest solution first.
    for v in range(v_max, 1, -1):
        if is_prime(v):
            u = 4 * v + 1
            if is_prime(u):
                p = 4 * u + 1
                # The check p < 1000000 is implicitly handled by starting v from v_max,
                # but we include it for clarity.
                if p < 1000000 and is_prime(p):
                    # Since we iterate downwards, the first solution found is the largest.
                    # The final output prints each number p, u, and v.
                    print(f"{p}:{u}:{v}")
                    return

    print("No solution found within the given constraints.")

if __name__ == '__main__':
    find_largest_prime_chain()