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

def find_largest_prime_chain():
    """
    Finds the largest prime p < 1,000,000 of the form p = 4u + 1,
    where u = 4v + 1, and p, u, v are all prime.
    """
    # The relationship is p = 4*(4*v + 1) + 1 = 16*v + 5.
    # To ensure p < 1,000,000, we have 16*v + 5 < 1,000,000.
    # 16*v < 999,995
    # v < 62499.6875
    # So, the maximum integer value for v is 62499.
    limit_v = 62499
    
    # We search downwards from the limit to find the largest p first.
    for v in range(limit_v, 1, -1):
        # Check if v is a prime number.
        if is_prime(v):
            u = 4 * v + 1
            # Check if u is a prime number.
            if is_prime(u):
                p = 4 * u + 1
                # Check if p is a prime number and within the machine's integer limit.
                if p < 1000000 and is_prime(p):
                    # Since we are iterating downwards, this is the largest triplet.
                    print(f"{p}:{u}:{v}")
                    # Print the equations as requested.
                    print(f"Equation for p: {p} = 4 * {u} + 1")
                    print(f"Equation for u: {u} = 4 * {v} + 1")
                    return

    print("No such prime triplet was found within the given constraints.")

# Run the function to find and print the result.
find_largest_prime_chain()