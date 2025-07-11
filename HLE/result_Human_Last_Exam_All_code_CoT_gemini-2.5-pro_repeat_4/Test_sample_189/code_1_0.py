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

def find_largest_prime_triplet():
    """
    Searches for the largest prime p such that p = 4u+1 and u = 4v+1,
    where p, u, and v are all prime, and p is less than 1,000,000.
    """
    # Initialize variables to store the largest found triplet
    largest_p = -1
    corresponding_u = -1
    corresponding_v = -1

    # The limit for v is derived from p = 16v + 5 < 1,000,000
    limit_v = 62500

    # Iterate through possible values for v
    for v in range(2, limit_v):
        # Check if v is a prime number
        if is_prime(v):
            # If v is prime, calculate u
            u = 4 * v + 1
            # Check if u is also prime
            if is_prime(u):
                # If u is prime, calculate p
                p = 4 * u + 1
                # Check if p is prime and within the machine's limits
                if p < 1000000 and is_prime(p):
                    # A valid triplet (p, u, v) is found.
                    # Since we are iterating v upwards, this will be the largest p found so far.
                    largest_p = p
                    corresponding_u = u
                    corresponding_v = v

    # After checking all candidates, print the result
    if largest_p != -1:
        # The problem asks to output each number in the final equation in the format p:u:v
        print(f"{largest_p}:{corresponding_u}:{corresponding_v}")
    else:
        print("No such prime triplet was found within the specified limits.")

# Execute the search and print the result
find_largest_prime_triplet()