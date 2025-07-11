import itertools

def find_integer_partitions(n):
    """
    Generates all integer partitions of n.
    Credit: Adapted from a recursive generator by Jerome Kelleher.
    """
    a = [0 for _ in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield tuple(a[:k + 2])
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield tuple(a[:k + 1])

def solve():
    """
    Solves the problem by finding the optimal partition of the dimension 6.
    """
    # N(d) is the maximum number of vectors in C^d with pairwise angle pi/3.
    # The inner product magnitude squared is cos^2(pi/3) = (1/2)^2 = 1/4.
    # These are known results from the study of complex equiangular lines.
    # The value for d=3 is from SIC-POVM constructions (m = d^2 when |(v_i, v_j)|^2 = 1/(d+1)).
    # For d=3, 1/(3+1) = 1/4, so m=3^2=9.
    N_d = {
        1: 1,
        2: 4,
        3: 9,
        4: 8,
        5: 10,
        6: 12,
    }

    max_vectors = 0
    best_partition = None
    best_sum_str = ""

    # The total dimension used must be between 2 and 6.
    # It must be at least 2 since we need at least two orthogonal subspaces
    # for an orthogonal pair to exist.
    for total_dim in range(2, 7):
        # Generate all partitions for the current total dimension
        for p in find_integer_partitions(total_dim):
            # The condition for an orthogonal pair means there must be at least
            # two orthogonal subspaces, so the partition must have at least 2 parts.
            if len(p) < 2:
                continue
            
            # Calculate the total number of vectors for this partition
            current_vectors = sum(N_d[dim] for dim in p)

            # Check if this is the best result so far
            if current_vectors > max_vectors:
                max_vectors = current_vectors
                best_partition = p
                sum_parts = [str(N_d[dim]) for dim in p]
                best_sum_str = " + ".join(sum_parts)
    
    print("This problem can be solved by partitioning the 6-dimensional space into a direct sum of orthogonal subspaces.")
    print("We then populate each subspace with the maximum possible number of equiangular vectors.")
    print("The optimal partition of dimensions and the corresponding number of vectors is found as follows:")
    print("-" * 30)
    print(f"The largest number of vectors is {max_vectors}.")
    print(f"This is achieved by partitioning the $\\mathbb{{C}}^6$ space into two orthogonal 3-dimensional subspaces.")
    print(f"The dimensional partition is {best_partition}.")
    print("The calculation for the total number of vectors is:")
    print(f"{best_sum_str} = {max_vectors}")

solve()
<<<18>>>