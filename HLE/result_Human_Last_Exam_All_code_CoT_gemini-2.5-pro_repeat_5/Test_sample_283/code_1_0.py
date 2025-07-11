import math

def solve_hss_multiplication_access():
    """
    Calculates the number of submatrices accessed during the multiplication
    of two matrices represented by a Hierarchical Semi-separable (HSS) tree.
    """
    # The depth of the HSS tree.
    d = 4

    # 1. Calculate the number of diagonal submatrices.
    # These correspond to the leaf nodes of the tree. At depth d, there are 2^d leaves.
    num_diag_blocks = 2**d

    # 2. Calculate the number of off-diagonal submatrices.
    # At each non-leaf level l (from 0 to d-1), there are 2^l parent nodes.
    # Each parent's children form a pair of siblings, corresponding to 2 off-diagonal blocks.
    off_diag_per_level = []
    for l in range(d):
        off_diag_per_level.append(2 * (2**l))
    num_off_diag_blocks = sum(off_diag_per_level)

    # 3. Calculate the total number of submatrices in a single HSS matrix representation.
    total_blocks_single_matrix = num_diag_blocks + num_off_diag_blocks

    # 4. Calculate the total number of accessed submatrices for multiplication (A * B).
    # This involves accessing the submatrices from both matrix A and matrix B.
    total_accessed_submatrices = 2 * total_blocks_single_matrix

    # --- Output the results step-by-step ---

    print(f"Analysis for an HSS tree of depth d = {d}:")
    print("-" * 40)

    # Print diagonal block calculation
    print("1. Calculation of Diagonal Submatrices (at leaf level):")
    print(f"   Number of diagonal submatrices = 2^{d} = {num_diag_blocks}")
    print()

    # Print off-diagonal block calculation
    print("2. Calculation of Off-Diagonal Submatrices (at non-leaf levels 0 to 3):")
    level_calcs = [f"(2 * 2^{l})" for l in range(d)]
    level_vals = [str(val) for val in off_diag_per_level]
    print(f"   Formula for each level l: 2 * 2^l")
    print(f"   Sum across levels = {' + '.join(level_calcs)}")
    print(f"                   = {' + '.join(level_vals)} = {num_off_diag_blocks}")
    print()

    # Print total for a single matrix
    print("3. Total submatrices in one HSS matrix:")
    print(f"   Total per matrix = (Diagonal) + (Off-Diagonal)")
    print(f"                    = {num_diag_blocks} + {num_off_diag_blocks} = {total_blocks_single_matrix}")
    print()

    # Print final calculation for multiplication
    print("4. Total submatrices accessed during multiplication (A * B):")
    print(f"   Total accessed = 2 * (Total per matrix)")
    print(f"                  = 2 * {total_blocks_single_matrix}")
    print(f"                  = {total_accessed_submatrices}")


solve_hss_multiplication_access()
