import math

def solve_hss_submatrices():
    """
    Calculates the number of submatrices accessed during matrix multiplication
    for a Hierarchical Semi-separable (HSS) tree of a given depth.

    The total number is the sum of:
    1. Off-diagonal blocks created at each level of partitioning.
    2. Diagonal blocks at the final leaf level.
    """
    # The depth of the Hierarchical Semi-separable tree.
    d = 4

    # --- 1. Calculate the number of off-diagonal blocks ---
    # The partitioning occurs at levels l = 0, 1, ..., d-1.
    # At each level l, there are 2^l nodes being partitioned, and each
    # partition creates 2 off-diagonal blocks.
    num_off_diagonals = 0
    off_diagonal_terms = []
    # Loop through the partitioning levels (0 to d-1)
    for l in range(d):
        term = 2**(l + 1)
        num_off_diagonals += term
        off_diagonal_terms.append(term)

    # --- 2. Calculate the number of leaf-level diagonal blocks ---
    # The recursion stops at the leaf level, d. The number of diagonal
    # blocks at this level is 2^d.
    num_leaf_diagonals = 2**d

    # --- 3. Calculate the total and print the equation ---
    total_submatrices = num_off_diagonals + num_leaf_diagonals

    # Construct the equation string to show each component
    off_diagonal_sum_str = " + ".join(map(str, off_diagonal_terms))
    
    print("A Hierarchical Semi-separable tree of depth 4 involves recursive partitioning.")
    print("We count the off-diagonal blocks from each partition level and the final diagonal blocks at the leaf level.")
    print("\nCalculation:")
    print(f"Off-diagonal blocks = (level 0) + (level 1) + (level 2) + (level 3)")
    print(f"                  = 2^(0+1) + 2^(1+1) + 2^(2+1) + 2^(3+1)")
    print(f"                  = {off_diagonal_terms[0]} + {off_diagonal_terms[1]} + {off_diagonal_terms[2]} + {off_diagonal_terms[3]} = {num_off_diagonals}")
    print(f"\nLeaf-level diagonal blocks = 2^4 = {num_leaf_diagonals}")
    print("\nTotal submatrices accessed = (Sum of off-diagonals) + (Leaf diagonals)")
    print("Final Equation:")
    print(f"{off_diagonal_sum_str} + {num_leaf_diagonals} = {total_submatrices}")


solve_hss_submatrices()
<<<46>>>