import collections

def count_allowed_pairs(limit):
    """
    Counts the number of allowed pairs (a,b) with 1 <= a, b <= limit.
    An ordered pair (a,b) is allowed if S(a) is a subset of S(b) or vice versa,
    where S(n) is the set of prime factors of n.
    """

    # Step 1: Compute prime factor sets for all numbers up to the limit.
    # We use a sieve-like method for efficiency.
    prime_factors = collections.defaultdict(set)
    for i in range(2, limit + 1):
        if not prime_factors[i]:  # i is prime
            for j in range(i, limit + 1, i):
                prime_factors[j].add(i)

    # Note: prime_factors[1] will be an empty set, which is correct as S(1) is empty.

    # Step 2: Initialize a counter
    allowed_pairs_count = 0

    # Step 3 & 4: Iterate through all pairs (a, b) and get their prime factor sets.
    for a in range(1, limit + 1):
        for b in range(1, limit + 1):
            s_a = prime_factors[a]
            s_b = prime_factors[b]

            # Step 5 & 6: Check the subset condition and increment the counter if it holds.
            if s_a.issubset(s_b) or s_b.issubset(s_a):
                allowed_pairs_count += 1
    
    # Step 7: Print the final result.
    print(f"The number of allowed pairs (a,b) with 1 <= a,b <= 1000 is:")
    print(allowed_pairs_count)

if __name__ == '__main__':
    # Set the limit for a and b
    N = 1000
    count_allowed_pairs(N)
