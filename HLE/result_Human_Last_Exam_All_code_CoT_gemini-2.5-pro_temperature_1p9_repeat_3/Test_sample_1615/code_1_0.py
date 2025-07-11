def solve_graph_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph G
    with 12345 vertices, given that G is not a complete graph.
    """
    
    # The number of vertices in the graph G
    n = 12345
    
    # The problem concerns the chromatic number, chi(G), of a graph G.
    # The graph G has n vertices and is not the complete graph K_n.

    # Step 1: Establish the upper bound for the chromatic number.
    # A standard theorem in graph theory states that any graph G with n vertices
    # that is not a complete graph can be colored with at most n-1 colors.
    # Proof idea: Since G is not complete, there exist two non-adjacent vertices.
    # These two vertices can be identified (or merged), resulting in a graph with n-1 vertices.
    # A graph with n-1 vertices can be colored with at most n-1 colors.
    # This coloring can be mapped back to the original graph G, showing chi(G) <= n-1.
    
    upper_bound = n - 1

    # Step 2: Show that this upper bound is achievable.
    # To do this, we need to construct a graph G' on n vertices that is not complete
    # and has a chromatic number of exactly n-1.
    
    # Consider the graph G' formed by the disjoint union of a complete graph on (n-1)
    # vertices (a clique K_{n-1}) and a single isolated vertex (K_1).
    # This graph has (n-1) + 1 = n vertices.
    # It is not a complete graph because the isolated vertex is not connected to any other vertex.
    
    # Step 3: Calculate the chromatic number of this example graph.
    # The chromatic number of a disjoint union of graphs is the maximum of the
    # chromatic numbers of its components.
    # chi(G') = max(chi(K_{n-1}), chi(K_1))
    
    # The chromatic number of a complete graph K_k is k.
    chi_K_n_minus_1 = n - 1
    
    # The chromatic number of a graph with one vertex K_1 is 1.
    chi_K_1 = 1
    
    # Calculate chi(G')
    chi_G_prime = max(chi_K_n_minus_1, chi_K_1)
    
    # Since we found a non-complete graph G' on n vertices with chi(G') = n-1,
    # the maximum possible chromatic number is indeed n-1.
    
    max_colors = n - 1

    print(f"The number of vertices is n = {n}.")
    print("The graph is not complete, so its chromatic number chi(G) is at most n-1.")
    print(f"To show this maximum is achievable, consider the graph K_({n}-1) + K_1.")
    print(f"This graph has a clique of size {n-1}, so it requires at least {n-1} colors.")
    print(f"Therefore, the maximum number of colors needed is {n} - 1 = {max_colors}.")

solve_graph_coloring_problem()