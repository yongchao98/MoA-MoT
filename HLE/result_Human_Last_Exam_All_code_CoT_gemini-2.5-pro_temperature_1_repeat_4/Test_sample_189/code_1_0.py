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

def find_special_primes():
    """
    Finds the largest prime p < 999999 of the form p=4u+1, where u=4v+1,
    and v, u, p are all prime numbers.
    """
    max_p_found = -1
    result_u = -1
    result_v = -1

    # From p = 4u+1 < 999999, we get u < 249999.5
    # From u = 4v+1 < 249999.5, we get v < 62499.625
    limit_v = 62500

    for v in range(2, limit_v):
        # Condition 1: v must be prime
        if is_prime(v):
            u = 4 * v + 1
            
            # Optimization: if u is already too large, no need to continue
            if u >= 250000:
                break

            # Condition 2: u must be prime
            if is_prime(u):
                p = 4 * u + 1

                # Condition 3: p must be prime and within the memory limit
                if p < 999999 and is_prime(p):
                    # A valid triplet is found. Check if it's the largest p.
                    if p > max_p_found:
                        max_p_found = p
                        result_u = u
                        result_v = v

    if max_p_found != -1:
        # Output the largest prime p found and the corresponding u and v
        print(f"{max_p_found}:{result_u}:{result_v}")
    else:
        print("No such prime satisfying the conditions was found.")

# Execute the search function
find_special_primes()