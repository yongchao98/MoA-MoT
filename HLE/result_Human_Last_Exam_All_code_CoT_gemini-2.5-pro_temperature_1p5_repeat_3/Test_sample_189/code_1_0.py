import math

def find_largest_prime_chain():
    """
    This function searches for the largest prime p < 1,000,000 such that
    p = 4u + 1 and u = 4v + 1, where p, u, and v are all prime numbers.
    """
    # The maximum value is determined by the ancient computer's architecture.
    LIMIT = 1000000

    # Use the Sieve of Eratosthenes to efficiently find all primes up to LIMIT.
    is_prime = [True] * LIMIT
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(LIMIT)) + 1):
        if is_prime[i]:
            for multiple in range(i * i, LIMIT, i):
                is_prime[multiple] = False

    # Initialize variables to store the largest found prime p and its components.
    max_p_found = -1
    result_u = -1
    result_v = -1

    # Determine the search limit for v.
    # If p < 1,000,000, then u < 250,000, and v < 62,500.
    v_limit = 62500

    # Iterate through all possible values for v to find the chain.
    for v in range(2, v_limit):
        # Check if v is a prime number.
        if is_prime[v]:
            u = 4 * v + 1
            # Check if u is within bounds and is a prime number.
            if u < LIMIT and is_prime[u]:
                p = 4 * u + 1
                # Check if p is within bounds and is a prime number.
                if p < LIMIT and is_prime[p]:
                    # A valid chain is found. Since we iterate v upwards,
                    # the last found p will be the largest.
                    max_p_found = p
                    result_u = u
                    result_v = v

    # Print the final result in the requested format.
    if max_p_found != -1:
        print(f"The largest prime found is p = {max_p_found}.")
        print(f"This is derived from the primes u = {result_u} and v = {result_v}.")
        print("The final equations are:")
        print(f"{max_p_found} = 4 * {result_u} + 1")
        print(f"{result_u} = 4 * {result_v} + 1")
    else:
        print("No such prime chain was found within the specified limits.")

# Execute the function to find and print the result.
find_largest_prime_chain()