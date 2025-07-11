def check_planarity_k33():
    """
    Checks if the K(3,3) graph (from the three utilities problem) can be planar.
    """
    # For the K(3,3) graph:
    # There are 3 houses and 3 utilities.
    num_houses = 3
    num_utilities = 3
    
    # The total number of vertices (V) is the sum of houses and utilities.
    V = num_houses + num_utilities
    
    # The number of edges (E) is the number of connections.
    # Each house connects to each utility.
    E = num_houses * num_utilities
    
    print(f"Analyzing the Three Utilities Problem as a graph K(3,3):")
    print(f"Number of vertices (V) = {V} (3 houses + 3 utilities)")
    print(f"Number of edges (E) = {E} (3 houses * 3 utilities)")
    print("-" * 30)
    
    # A key property of planar graphs with no triangles (cycles of length 3) is E <= 2V - 4.
    # K(3,3) is a bipartite graph, so it has no odd-length cycles, and thus no triangles.
    # We can check if this condition holds.
    print("For a graph with no triangles to be planar, it must satisfy the inequality: E <= 2*V - 4")
    
    # Perform the check
    condition_holds = E <= 2 * V - 4
    
    # Print the equation with the actual numbers
    print("\nSubstituting the values into the formula:")
    print(f"{E} <= 2 * {V} - 4")
    
    # Simplify the right side of the equation
    rhs = 2 * V - 4
    print(f"{E} <= {rhs}")
    
    print("\nResult:")
    if condition_holds:
        print("The condition is met. The graph could be planar.")
    else:
        print("The statement is FALSE. The condition for planarity is not met.")
        print("This proves that the K(3,3) graph is non-planar.")

check_planarity_k33()
<<<E>>>