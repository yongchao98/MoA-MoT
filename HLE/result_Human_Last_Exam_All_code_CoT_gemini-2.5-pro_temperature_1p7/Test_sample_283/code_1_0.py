import math

def calculate_hss_submatrices(depth):
    """
    Calculates the number of submatrices in an HSS representation of a given depth.
    
    Args:
        depth (int): The depth of the HSS tree.
    """
    print(f"Analyzing an HSS tree of depth d = {depth}.")
    print("The total number of submatrices is the sum of generator matrices for off-diagonal blocks and dense matrices for leaf-level diagonal blocks.")
    print("-" * 50)
    
    # --- Step 1: Calculate the number of generator matrices ---
    print("Step 1: Counting generator matrices from off-diagonal blocks.")
    print("At each partitioning level k (from 1 to d), 2^k off-diagonal blocks are formed.")
    print("Each off-diagonal block is stored as two 'generator' matrices (U and V).")
    
    num_generators = 0
    generator_terms_str = []
    
    for k in range(1, depth + 1):
        off_diagonal_blocks = 2**k
        generators_at_level_k = off_diagonal_blocks * 2
        num_generators += generators_at_level_k
        generator_terms_str.append(str(generators_at_level_k))
        print(f"  - Level {k}: 2^{k} = {off_diagonal_blocks} off-diagonal blocks => {generators_at_level_k} generator matrices.")

    print(f"Total generator matrices = {' + '.join(generator_terms_str)} = {num_generators}")
    print("-" * 50)
    
    # --- Step 2: Calculate the number of leaf blocks ---
    print("Step 2: Counting the diagonal blocks at the leaf level.")
    print("After d levels of partitioning, there are 2^d final diagonal blocks (leaves).")
    
    num_leaves = 2**depth
    print(f"Number of leaf blocks = 2^{depth} = {num_leaves}")
    print("-" * 50)

    # --- Step 3: Sum everything up ---
    print("Step 3: Calculating the total number of submatrices.")
    total_submatrices = num_generators + num_leaves
    
    # Display the final equation with all numeric components.
    final_equation = f"Total submatrices = ({' + '.join(generator_terms_str)}) + {num_leaves} = {total_submatrices}"
    
    print("The final calculation is:")
    print(final_equation)


# Set the depth for the problem
tree_depth = 4
calculate_hss_submatrices(tree_depth)