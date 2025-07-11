def solve_vector_problem():
    """
    Solves the problem by finding the optimal partition of the 6-dimensional space.
    """
    
    # k_map stores the maximum number of vectors k for a given subspace dimension d.
    # These are based on known results for complex equiangular lines with |(v,w)|=1/2.
    k_map = {
        1: 1,
        2: 3,
        3: 9
    }

    # The total dimension of the vector space.
    total_dimension = 6
    
    # We are looking for partitions of 6 into integers from our k_map keys {1, 2, 3}.
    # The problem requires at least one orthogonal pair, so we need at least two subspaces in our partition.
    partitions = [
        [3, 3],
        [3, 2, 1],
        [3, 1, 1, 1],
        [2, 2, 2],
        [2, 2, 1, 1],
        [2, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1]
    ]

    max_vectors = 0
    best_partition = []
    
    print("Evaluating possible partitions of the 6-dimensional space:")
    print("-" * 60)

    for p in partitions:
        current_vectors = 0
        
        # Build the equation string for printing
        dim_str = " + ".join(map(str, p))
        k_values = [k_map[d] for d in p]
        k_str = " + ".join(map(str, k_values))
        current_vectors = sum(k_values)

        print(f"Partition {dim_str} = {total_dimension}:")
        print(f"Number of vectors = {k_str} = {current_vectors}")
        
        if current_vectors > max_vectors:
            max_vectors = current_vectors
            best_partition = p
        print("-" * 20)

    print("\nSummary:")
    best_dim_str = " + ".join(map(str, best_partition))
    best_k_values = [k_map[d] for d in best_partition]
    best_k_str = " + ".join(map(str, best_k_values))

    print(f"The optimal strategy is to partition the C^6 space into two orthogonal 3-dimensional subspaces ({best_dim_str}).")
    print(f"This configuration allows for a total of {best_k_str} = {max_vectors} vectors.")
    print(f"The largest number of such vectors is {max_vectors}.")

solve_vector_problem()