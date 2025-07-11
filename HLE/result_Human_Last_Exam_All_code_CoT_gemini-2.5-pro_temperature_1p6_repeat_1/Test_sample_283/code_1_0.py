import math

def solve_hss_submatrices():
    """
    Calculates the number of submatrices in an HSS tree of a given depth.

    In a Hierarchical Semi-separable (HSS) representation, a matrix is
    recursively partitioned. This structure is represented by a binary tree.
    - At each internal node, the matrix block is partitioned into 2x2 sub-blocks.
      The two off-diagonal blocks are compressed, creating 2 submatrix representations.
    - At each leaf node, the corresponding diagonal block is stored as a
      single dense submatrix.

    The total number of submatrices is the sum of all compressed off-diagonal
    blocks from all internal nodes and all dense diagonal blocks from leaf nodes.
    """
    # The depth of the HSS tree. The root is at level 0, leaves are at level 4.
    depth = 4

    # Internal nodes are at levels 0, 1, ..., depth-1.
    # The number of internal nodes is the sum of a geometric series: 2^0 + 2^1 + ... + 2^(depth-1).
    # This sum equals 2^depth - 1.
    num_internal_nodes = int(math.pow(2, depth)) - 1

    # Each internal node defines 2 compressed off-diagonal submatrices.
    num_off_diagonal_blocks = 2 * num_internal_nodes

    # Leaf nodes are at level `depth`.
    # The number of leaf nodes is 2^depth.
    num_leaf_nodes = int(math.pow(2, depth))

    # Each leaf node corresponds to one dense diagonal submatrix.
    num_diagonal_leaf_blocks = num_leaf_nodes

    # The total number of submatrices accessed during a full matrix operation
    # is the sum of these components.
    total_submatrices = num_off_diagonal_blocks + num_diagonal_leaf_blocks

    print(f"For an HSS tree of depth {depth}:")
    print(f"1. The number of internal nodes is 2^{depth} - 1 = {num_internal_nodes}.")
    print(f"   Each represents 2 off-diagonal blocks, so there are {num_off_diagonal_blocks} off-diagonal submatrices.")
    print(f"2. The number of leaf nodes is 2^{depth} = {num_leaf_nodes}.")
    print(f"   Each represents 1 diagonal block, so there are {num_diagonal_leaf_blocks} diagonal submatrices.")
    print("\nThe total number of submatrices is the sum of these two counts.")
    print(f"Total Submatrices = {num_off_diagonal_blocks} + {num_diagonal_leaf_blocks} = {total_submatrices}")

solve_hss_submatrices()