import math

def is_prime(n):
    """
    Checks if a number n is prime using an efficient algorithm.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    # Check for divisors of the form 6k +/- 1 up to sqrt(n)
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_special_primes():
    """
    Finds the largest prime p < 999999 such that p = 4u+1 and u = 4v+1,
    where u and v are also primes.
    """
    largest_p_found = -1
    corresponding_u = -1
    corresponding_v = -1

    # The maximum value for p is 999999.
    # p = 16v + 5 < 999999 => 16v < 999994 => v < 62499.625
    v_limit = 62499

    for v in range(2, v_limit + 1):
        # 1. Check if v is prime
        if is_prime(v):
            # 2. Calculate u and check if u is prime
            u = 4 * v + 1
            if is_prime(u):
                # 3. Calculate p and check if p is prime
                p = 4 * u + 1
                # The p < 999999 check is implicitly handled by v_limit,
                # but we keep it for clarity.
                if p < 999999 and is_prime(p):
                    # Since we are iterating upwards, the last one found will be the largest.
                    largest_p_found = p
                    corresponding_u = u
                    corresponding_v = v

    if largest_p_found != -1:
        # Output the final result in the format p:u:v
        print(f"{largest_p_found}:{corresponding_u}:{corresponding_v}")
    else:
        print("No such prime triple found within the given limits.")

if __name__ == "__main__":
    find_special_primes()