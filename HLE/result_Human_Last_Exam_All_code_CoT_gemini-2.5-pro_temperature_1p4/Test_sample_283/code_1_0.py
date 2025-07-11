import math

def solve_hss_submatrices():
    """
    Calculates the number of submatrices accessed during a matrix multiplication
    for a Hierarchical Semi-separable (HSS) tree of a given depth.
    """
    # The depth of the Hierarchical Semi-separable tree.
    depth = 4

    # 1. Calculate the number of diagonal submatrices at the leaf level.
    # For a binary tree of depth 'd' (levels 0 to d-1), there are 2^(d-1) leaves.
    num_diagonal_blocks = 2**(depth - 1)

    # 2. Calculate the number of off-diagonal submatrices.
    # There are 2 off-diagonal blocks for each non-leaf node. The number of
    # non-leaf nodes in a complete binary tree of depth 'd' is 2^(d-1) - 1.
    num_non_leaf_nodes = 2**(depth - 1) - 1
    num_off_diagonal_blocks = 2 * num_non_leaf_nodes

    # 3. Calculate the total number of submatrices accessed.
    total_submatrices = num_diagonal_blocks + num_off_diagonal_blocks

    # 4. Print the results in a clear equation format.
    print(f"For an HSS tree of depth d = {depth}:")
    print(f"Number of diagonal leaf submatrices accessed = 2^(d-1) = 2^({depth-1}) = {num_diagonal_blocks}")
    print(f"Number of off-diagonal submatrices accessed = 2 * (2^(d-1) - 1) = 2 * (2^({depth-1}) - 1) = {num_off_diagonal_blocks}")
    print("\nFinal Equation:")
    print(f"Total submatrices = {num_diagonal_blocks} (diagonal) + {num_off_diagonal_blocks} (off-diagonal) = {total_submatrices}")

solve_hss_submatrices()