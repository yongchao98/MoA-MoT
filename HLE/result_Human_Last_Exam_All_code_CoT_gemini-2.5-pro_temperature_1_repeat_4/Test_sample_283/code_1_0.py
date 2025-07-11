def solve_hss_multiplication():
    """
    Calculates the number of submatrices accessed during a matrix multiplication
    for a Hierarchical Semi-separable (HSS) tree of a given depth.
    """
    depth = 4

    # An HSS tree of depth D has D+1 levels, from 0 to D.
    levels = list(range(depth + 1))

    # The number of nodes (submatrices) at level l is 2^l.
    nodes_per_level = [2**l for l in levels]

    # The total number of submatrices in one HSS tree is the sum of nodes at all levels.
    nodes_per_tree = sum(nodes_per_level)

    # Matrix multiplication involves two matrices (A and B).
    total_submatrices_accessed = nodes_per_tree * 2

    # Print the explanation and the final equation.
    print(f"A Hierarchical Semi-separable (HSS) tree of depth {depth} has levels from 0 to {depth}.")
    print(f"The number of submatrices at each level is: {nodes_per_level}")

    # Create the string for the sum equation for one tree
    sum_str = " + ".join(map(str, nodes_per_level))
    print(f"The total number of submatrices in a single HSS tree is the sum of nodes at all levels:")
    print(f"{sum_str} = {nodes_per_tree}")
    
    print("\nMatrix multiplication involves two such trees (one for each matrix).")
    print("Therefore, the total number of submatrices accessed is:")
    print(f"{nodes_per_tree} (from matrix A) + {nodes_per_tree} (from matrix B) = {total_submatrices_accessed}")

solve_hss_multiplication()
