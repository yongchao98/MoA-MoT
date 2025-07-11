def solve_reversal_puzzle():
    """
    Calculates the minimum moves to reverse a sequence of 100 elements
    with special swap rules.
    """
    N = 100  # Total number of elements
    k = 5    # The non-adjacent swap distance is i -> i+5

    # The problem reduces to sorting k "buckets" of elements.
    # Each bucket contains N/k elements.
    bucket_size = N // k
    print(f"The {N} positions are divided into {k} buckets based on their index modulo {k}.")
    print(f"Each bucket contains {N}/{k} = {bucket_size} elements.\n")

    # Determine the target permutation of the buckets.
    # An element at position i moves to 101-i.
    # A bucket r (i % k) moves to bucket (101-i) % k = (1-r) % k.
    # Let's create this permutation for r = 0, 1, ..., k-1.
    # Using 1-based indexing for explanation, 0-based for code.
    perm = [(1 - r + k) % k for r in range(k)]
    
    print(f"The initial order of populations in buckets [0, 1, 2, 3, 4] is (S0, S1, S2, S3, S4).")
    print(f"The target order is (S{perm[0]}, S{perm[1]}, S{perm[2]}, S{perm[3]}, S{perm[4]}).")
    print(f"This corresponds to sorting the permutation: {perm}\n")

    # The number of moves is bucket_size * (number of adjacent swaps to sort the permutation).
    # The number of adjacent swaps to sort a permutation is its number of inversions.
    inversions = 0
    inversion_pairs = []
    for i in range(k):
        for j in range(i + 1, k):
            if perm[i] > perm[j]:
                inversions += 1
                inversion_pairs.append((perm[i], perm[j]))

    print(f"The number of inversions in the permutation {perm} is {inversions}.")
    # print(f"The inversions are: {inversion_pairs}")
    
    # Each adjacent swap of buckets costs `bucket_size` moves.
    total_moves = inversions * bucket_size
    
    print(f"\nFinal Calculation:")
    print(f"Minimum Moves = (Number of Inversions) * (Elements per Bucket)")
    print(f"Minimum Moves = {inversions} * {bucket_size}")
    print(f"Minimum Moves = {total_moves}")

solve_reversal_puzzle()
<<<80>>>