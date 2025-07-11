def solve_hss_multiplication():
    """
    Calculates the number of submatrices accessed during the multiplication
    of two HSS matrices represented by a tree of a given depth.
    """
    # The depth of the Hierarchical Semi-separable (HSS) tree.
    depth = 4

    # 1. Calculate the number of nodes in the binary tree.
    # A full binary tree of depth 'd' has 2**(d+1) - 1 total nodes.
    total_nodes = 2**(depth + 1) - 1
    
    # 2. Calculate the number of leaf nodes.
    # These correspond to the dense diagonal 'D' blocks.
    num_leaf_nodes = 2**depth
    num_d_blocks = num_leaf_nodes

    # 3. Calculate the number of 'U' and 'V' generator blocks.
    # These are defined for all nodes in the tree except for the root node.
    num_uv_nodes = total_nodes - 1
    num_u_blocks = num_uv_nodes
    num_v_blocks = num_uv_nodes

    # 4. Calculate the number of 'B' translation blocks.
    # These are defined for each pair of siblings, which is equivalent to
    # the number of non-leaf nodes in the tree.
    num_non_leaf_nodes = total_nodes - num_leaf_nodes
    num_b_blocks = num_non_leaf_nodes
    
    # 5. Calculate the total number of submatrices for a single HSS matrix.
    total_per_matrix = num_d_blocks + num_u_blocks + num_v_blocks + num_b_blocks

    # 6. For a matrix multiplication (A * B), the submatrices of both matrices are accessed.
    total_accessed = 2 * total_per_matrix

    # --- Output the results step-by-step ---
    print(f"Calculation for an HSS matrix defined by a tree of depth d = {depth}:")
    print("-" * 60)
    
    # Print count for D blocks
    print(f"Number of 'D' blocks (one for each leaf node):")
    print(f"  2^d = 2^{depth} = {num_d_blocks}")
    print()

    # Print count for U and V blocks
    print(f"Number of 'U' and 'V' blocks (one for each non-root node):")
    print(f"  Total nodes = 2^(d+1) - 1 = {total_nodes}")
    print(f"  Non-root nodes = {total_nodes} - 1 = {num_uv_nodes}")
    print(f"  Number of 'U' blocks = {num_u_blocks}")
    print(f"  Number of 'V' blocks = {num_v_blocks}")
    print()

    # Print count for B blocks
    print(f"Number of 'B' blocks (one for each non-leaf node):")
    print(f"  Non-leaf nodes = Total nodes - Leaf nodes = {total_nodes} - {num_leaf_nodes} = {num_non_leaf_nodes}")
    print(f"  Number of 'B' blocks = {num_b_blocks}")
    print()

    print("-" * 60)
    print("Total submatrices for a single HSS matrix:")
    print(f"  N_one = (D blocks) + (U blocks) + (V blocks) + (B blocks)")
    print(f"  N_one = {num_d_blocks} + {num_u_blocks} + {num_v_blocks} + {num_b_blocks} = {total_per_matrix}")
    print()
    
    print("Total submatrices accessed during matrix multiplication (A * B):")
    print(f"  N_total = 2 * N_one")
    print(f"  The final equation is:")
    print(f"  Total = 2 * ({num_d_blocks} + {num_u_blocks} + {num_v_blocks} + {num_b_blocks}) = {total_accessed}")

solve_hss_multiplication()