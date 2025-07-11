def solve_graph_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a non-complete
    graph with a given number of vertices.
    """
    
    # The number of vertices in the graph G
    num_vertices = 12345

    # Step-by-step explanation
    print("This problem asks for the maximum chromatic number (χ(G)) of a graph G with 12345 vertices, given that G is not a complete graph.")
    print("\nHere's the reasoning:")
    print("1. For any graph with n vertices, the maximum possible chromatic number is n. This value is achieved only by the complete graph, K_n, where every vertex is connected to every other.")
    print(f"2. Our graph G has n = {num_vertices} vertices but is not a complete graph.")
    print("3. A theorem in graph theory states that for any graph with n vertices that is not complete, its chromatic number is at most n-1.")
    print("4. This maximum value of n-1 is achievable. For example, a graph consisting of a complete subgraph on n-1 vertices (a K_{n-1} clique) and one additional isolated vertex has a chromatic number of exactly n-1.")
    print("\nTherefore, the maximum number of colors needed is n - 1.")

    # Calculate the result
    max_colors = num_vertices - 1

    # Print the final calculation as an equation
    print("\nFinal Calculation:")
    print(f"{num_vertices} - 1 = {max_colors}")

solve_graph_coloring_problem()