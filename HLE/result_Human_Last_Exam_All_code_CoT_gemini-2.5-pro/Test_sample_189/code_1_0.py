def solve_prime_chain():
    """
    Finds the largest prime p < 1,000,000 of the form p = 4u + 1,
    where u is a prime of the form u = 4v + 1, and v is also prime.
    """
    limit = 1000000

    # Step 1: Use a Sieve of Eratosthenes for efficient prime checking.
    is_prime = [True] * limit
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for multiple in range(i * i, limit, i):
                is_prime[multiple] = False

    max_p = 0
    final_u = 0
    final_v = 0

    # Step 2: Determine the search range for v.
    # If p = 4u+1 < 1M -> u < 250k.
    # If u = 4v+1 < 250k -> v < 62.5k.
    limit_v = 62500

    # Step 3: Iterate through possible values for v and check the chain.
    for v_candidate in range(2, limit_v):
        if is_prime[v_candidate]:
            # We have a prime v. Calculate u.
            u_candidate = 4 * v_candidate + 1

            if u_candidate < limit and is_prime[u_candidate]:
                # We have a prime u of the correct form. Calculate p.
                p_candidate = 4 * u_candidate + 1

                if p_candidate < limit and is_prime[p_candidate]:
                    # We found a valid (p, u, v) triplet.
                    # Check if this p is the largest found so far.
                    if p_candidate > max_p:
                        max_p = p_candidate
                        final_u = u_candidate
                        final_v = v_candidate

    # Step 4: Print the final results.
    if max_p > 0:
        print(f"The largest prime p found is: {max_p}")
        print("It was generated from the following prime chain:")
        print(f"{final_u} = 4 * {final_v} + 1")
        print(f"{max_p} = 4 * {final_u} + 1")
        print("\nFinal Answer:")
        print(f"{max_p}:{final_u}:{final_v}")
    else:
        print("No prime chain p = 4u+1 (u=4v+1) was found within the limit of 1,000,000.")

solve_prime_chain()