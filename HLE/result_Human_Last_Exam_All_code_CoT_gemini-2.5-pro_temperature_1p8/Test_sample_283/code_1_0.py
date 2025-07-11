import math

def solve_hss_submatrices():
    """
    Calculates the number of submatrices accessed during matrix multiplication
    for an HSS tree of a given depth.
    """
    # The depth of the Hierarchical Semi-separable tree.
    depth = 4

    print("For a Hierarchical Semi-separable (HSS) tree, matrix multiplication involves accessing data associated with each node.")
    print(f"The tree has a depth of {depth}, meaning it has levels from 0 to {depth}.")
    print("The number of nodes at each level 'k' in this binary tree is 2^k.")
    print("To find the total number of submatrices accessed, we sum the nodes at each level.")
    print("\nCalculation:")
    
    # Calculate the number of nodes at each level from 0 to depth.
    nodes_per_level = []
    for k in range(depth + 1):
        nodes_per_level.append(int(math.pow(2, k)))
        
    # Calculate the total number of nodes.
    total_nodes = sum(nodes_per_level)
    
    # Create the equation string to display the calculation.
    # e.g., "1 + 2 + 4 + ..."
    equation_str = " + ".join(map(str, nodes_per_level))
    
    print(f"Total Submatrices = {equation_str} = {total_nodes}")

solve_hss_submatrices()