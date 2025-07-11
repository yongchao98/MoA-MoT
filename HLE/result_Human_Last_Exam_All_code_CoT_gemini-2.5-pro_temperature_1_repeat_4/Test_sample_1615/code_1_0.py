def solve_graph_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph
    with 12345 vertices that is not a complete graph.
    """
    n = 12345

    print("Let G be a graph with n vertices, where n = 12345.")
    print("We are given that G is not a complete graph (G != K_n).")
    print("The question asks for the maximum number of colors needed for a proper vertex coloring of G.")
    print("This is equivalent to finding the maximum possible value of the chromatic number, chi(G), for such a graph G.")

    print("\n--- Step 1: Finding the upper bound ---")
    print("Since G is not a complete graph, there must be at least one pair of vertices, let's call them 'u' and 'v', that are not connected by an edge.")
    print("Because u and v are not adjacent, they can be assigned the same color in a proper coloring.")
    print("Let's use Color 1 for both u and v.")
    print(f"The remaining n - 2 = {n - 2} vertices can be colored using new, distinct colors.")
    print(f"This coloring strategy uses a total of 1 (for u and v) + {n - 2} (for the rest) = {n - 1} colors.")
    print(f"This shows that any non-complete graph with {n} vertices can be colored with at most {n - 1} colors.")
    print(f"Therefore, the maximum chromatic number is at most {n - 1}.")

    print("\n--- Step 2: Showing the upper bound is achievable ---")
    print("To show that this bound can be reached, we need to find an example of a non-complete graph on n vertices that requires exactly n-1 colors.")
    print(f"Consider the graph G* created by taking a complete graph on {n} vertices (K_{n}) and removing a single edge.")
    print("This graph G* is not complete, so it meets the problem's conditions.")
    print(f"Let the removed edge be between vertices 'u' and 'v'.")
    print(f"Now, consider the subgraph of G* induced by the set of all vertices except 'v'. This subgraph has {n-1} vertices.")
    print("In this subgraph, all vertices are mutually adjacent, forming a complete graph K_{n-1} (also known as an (n-1)-clique).")
    print(f"A complete graph on {n-1} vertices, K_{n-1}, requires exactly {n-1} colors for a proper coloring.")
    print("Since G* contains a subgraph that requires {n-1} colors, G* itself must require at least {n-1} colors.")
    
    print("\n--- Conclusion ---")
    result = n - 1
    print(f"From Step 1, we know the maximum number of colors is at most {n - 1}.")
    print(f"From Step 2, we found a graph that requires at least {n - 1} colors.")
    print(f"Combining these, the maximum number of colors needed is exactly {n-1}.")
    print("\nFor n = 12345, the final calculation is:")
    print(f"{n} - 1 = {result}")

solve_graph_coloring_problem()