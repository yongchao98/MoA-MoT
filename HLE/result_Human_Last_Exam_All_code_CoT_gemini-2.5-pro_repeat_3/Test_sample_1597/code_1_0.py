def solve_allowed_pairs():
    """
    Calculates the number of allowed pairs (a, b) in a given range.
    A pair (a, b) is allowed if the set of prime factors of a is a subset of
    the set of prime factors of b, or vice-versa.
    """
    limit = 1000

    # Step 1: Pre-compute the set of prime factors for each number up to the limit.
    # We use a sieve-like method. For each prime we find, we add it to the factor
    # sets of all its multiples.
    prime_factor_sets = [set() for _ in range(limit + 1)]
    for i in range(2, limit + 1):
        # If i has no prime factors marked yet, it must be a prime number.
        if not prime_factor_sets[i]:
            # i is a prime. Mark all its multiples.
            for j in range(i, limit + 1, i):
                prime_factor_sets[j].add(i)

    # Step 2: Iterate through all pairs and count the allowed ones.
    allowed_pairs_count = 0
    for a in range(1, limit + 1):
        for b in range(1, limit + 1):
            # Step 3: Retrieve the pre-computed sets of prime factors.
            set_a = prime_factor_sets[a]
            set_b = prime_factor_sets[b]

            # Step 4: Check if one set is a subset of the other.
            if set_a.issubset(set_b) or set_b.issubset(set_a):
                allowed_pairs_count += 1

    # Step 5: Print the final count.
    print(allowed_pairs_count)

solve_allowed_pairs()