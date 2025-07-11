def solve_coloring_problem():
    """
    This function calculates the maximum number of colors needed for a proper vertex
    coloring of a graph with 12345 vertices that is not a complete graph.
    """
    
    # Number of vertices in the graph G
    n = 12345

    # Explain the reasoning behind the solution.
    print("This problem asks for the maximum number of colors needed to properly color a graph G with 12345 vertices, given that G is not a complete graph.")
    print("Let's break down the logic:")
    
    print("\nStep 1: The absolute maximum chromatic number")
    print(f"For any graph with n = {n} vertices, the maximum possible number of colors for a proper coloring is n.")
    print("This maximum is only achieved if the graph is a complete graph (K_n), where every vertex is connected to every other vertex.")

    print("\nStep 2: The constraint of the graph not being complete")
    print("The problem states that G is not a complete graph. This means there is at least one pair of vertices that are NOT connected by an edge.")
    print("Since these two vertices are not adjacent, they can be assigned the same color. This fact implies that the chromatic number of G must be less than n.")
    print(f"Thus, for any non-complete graph G with {n} vertices, its chromatic number chi(G) is at most n - 1.")

    print("\nStep 3: Finding a graph that meets this maximum")
    print("To confirm that n-1 is the maximum possible value, we must show that it is achievable.")
    print(f"Consider a graph constructed by taking the complete graph on {n} vertices, K_{n}, and removing a single edge. This graph is not complete.")
    print(f"This graph still contains a clique (a fully connected subgraph) of size n-1. A clique of size n-1 requires n-1 distinct colors.")
    print(f"Therefore, there exists a non-complete graph with {n} vertices that requires n-1 colors.")

    print("\nStep 4: Conclusion and Calculation")
    print("From the logic above, the maximum number of colors needed for a non-complete graph on n vertices is exactly n-1.")
    
    # Perform the final calculation
    max_colors = n - 1
    
    print("\nPerforming the calculation for n = 12345:")
    print("The final equation is:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()