import math

def solve_hss_submatrices():
    """
    Calculates the number of submatrices in an HSS tree of a given depth.
    """
    # The depth of the Hierarchical Semi-separable tree.
    depth = 4

    print("For a Hierarchical Semi-separable tree, the total number of accessed submatrices")
    print("corresponds to the total number of nodes in the tree representation.")
    print(f"Given a tree of depth {depth}, we have levels 0 through {depth}.")
    print("The number of nodes (submatrices) at each level 'l' is 2^l.")
    print("\nTo find the total, we sum the number of nodes at each level:")

    # Calculate the number of nodes at each level
    nodes_per_level = []
    for level in range(depth + 1):
        num_nodes = int(math.pow(2, level))
        nodes_per_level.append(num_nodes)

    # Format the equation string
    equation_parts = [str(n) for n in nodes_per_level]
    equation_str = " + ".join(equation_parts)

    # Calculate the total sum
    total_nodes = sum(nodes_per_level)

    # Print the final equation and the result
    print(f"Total Submatrices = {equation_str}")
    print(f"Total Submatrices = {total_nodes}")

solve_hss_submatrices()