def find_prime_triplet():
    """
    This function finds the largest prime p such that p = 4u+1 and u = 4v+1,
    where u is also prime, for p up to a defined LIMIT.
    """
    # Set the search limit based on the problem's context.
    LIMIT = 1000000

    # Step 1: Generate primes using Sieve of Eratosthenes for efficient lookup.
    is_prime = [True] * (LIMIT + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(LIMIT**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i * i, LIMIT + 1, i):
                is_prime[multiple] = False

    # Create a list of all prime numbers found by the sieve.
    primes_list = [i for i, is_p in enumerate(is_prime) if is_p]

    # Step 2 & 3: Iterate backwards through the list of primes to find the largest p first.
    found = False
    for p in reversed(primes_list):
        # Condition 1: Check if p is of the form 4u + 1.
        # This is equivalent to checking if p = 1 (mod 4).
        if (p - 1) % 4 == 0:
            u = (p - 1) // 4

            # Condition 2: Check if the calculated u is also a prime number.
            # We must also ensure u is within the bounds of our sieve array.
            if u < len(is_prime) and is_prime[u]:
                
                # Condition 3: Check if prime u is of the form 4v + 1.
                # This is equivalent to checking if u = 1 (mod 4).
                if (u - 1) % 4 == 0:
                    v = (u - 1) // 4
                    
                    # Since we are iterating from the largest p, the first triplet
                    # found is the solution.
                    # The final print statement fulfills the request to output each number
                    # in the final equation via the p:u:v format.
                    print(f"{p}:{u}:{v}")
                    found = True
                    break

    if not found:
        print("No such prime triplet was found within the specified limit.")

# Execute the function to find and print the result.
find_prime_triplet()