def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph G with 12345 vertices, given that G is not a complete graph.
    """
    
    n = 12345
    
    print("Let n be the number of vertices in the graph G.")
    print(f"n = {n}")
    print("\nThe problem asks for the maximum number of colors needed for a proper coloring of G, which is the maximum possible value of the chromatic number, chi(G).")
    print("The graph G has n=12345 vertices and is not a complete graph (K_n).\n")
    
    print("Step 1: Find an upper bound for chi(G).")
    print("A complete graph K_n requires n colors, as every vertex is connected to every other vertex.")
    print("However, the problem states that G is NOT a complete graph.")
    print("This means there is at least one pair of vertices, let's call them u and v, that are not connected by an edge.")
    print("Because u and v are not adjacent, they can be assigned the same color in a proper coloring.")
    print("We can devise a coloring scheme: color u and v with color 1. Then, color the remaining n-2 vertices with n-2 new, distinct colors (color 2, color 3, ..., color n-1).")
    print("This is a valid proper coloring, and it uses (n-2) + 1 = n-1 colors.")
    print(f"This shows that for any non-complete graph G with n vertices, chi(G) <= n - 1 = {n-1}.\n")
    
    print("Step 2: Show that this upper bound is achievable.")
    print("To show that n-1 is the maximum possible number of colors, we need to show that there exists a graph G (which is not K_n) that requires exactly n-1 colors.")
    print("Consider a graph G' constructed by taking the complete graph K_n and removing a single edge, say between vertices u and v.")
    print("This graph G' is not complete, so it satisfies the condition of the problem.")
    print("Now, let's find the chromatic number of G'. Consider the set of vertices V' = V - {v}, where V is the set of all n vertices.")
    print("The subgraph induced by the n-1 vertices in V' contains a clique of size n-1 (a K_{n-1}) because u is connected to all vertices in V' except itself, and all other vertices in V' are mutually connected.")
    print("Wait, let's correct that argument for the clique.")
    print("Let the vertices be v1, v2, ..., vn. Let the removed edge be (v_{n-1}, v_n).")
    print("Consider the set of n-1 vertices {v1, v2, ..., v_{n-1}}.")
    print("In this set, every pair of vertices is connected by an edge in our graph G'. This forms a clique of size n-1.")
    print("Since G' contains a clique of size n-1, any proper coloring of G' requires at least n-1 colors. So, chi(G') >= n-1.\n")
    
    print("Step 3: Conclusion.")
    print("From Step 1, we know that chi(G) <= n-1.")
    print("From Step 2, we found a specific graph G' for which chi(G') >= n-1.")
    print("Therefore, the maximum possible chromatic number for a graph G with n vertices that is not a complete graph is exactly n-1.")
    
    result = n - 1
    
    print("\nFinal Calculation:")
    print(f"The maximum number of colours is n - 1.")
    print(f"Maximum colours = {n} - 1 = {result}")

solve_coloring_problem()