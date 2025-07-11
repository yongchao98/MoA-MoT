def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph with
    12345 vertices that is not a complete graph.
    """
    
    # Number of vertices in the graph G
    n = 12345
    
    print("Step 1: Understand the problem.")
    print(f"We have a graph G with n = {n} vertices.")
    print("The graph G is not a complete graph (K_n), where every vertex is connected to every other vertex.")
    print("We want to find the maximum possible number of colors needed for a proper vertex coloring, which is the maximum possible value of the chromatic number, chi(G).\n")
    
    print("Step 2: Find an upper bound for the number of colors.")
    print(f"For any graph with n vertices, the chromatic number chi(G) is at most n.")
    print(f"A graph needs n colors (i.e., chi(G) = n) if and only if it is the complete graph K_n.")
    print(f"Since G is not a complete graph, its chromatic number must be less than n.")
    print(f"So, chi(G) <= n - 1. For our case, chi(G) <= {n} - 1 = {n-1}.\n")
    
    print("Step 3: Show that this maximum value is achievable.")
    print("We need to check if a graph exists that satisfies the conditions and requires exactly n - 1 colors.")
    print("Consider a graph G' made of two disconnected components:")
    print(f"1. A complete graph on n - 1 = {n-1} vertices (K_{n-1}).")
    print("2. A single isolated vertex (K_1).")
    print(f"This graph G' has a total of ({n-1}) + 1 = {n} vertices.")
    print("G' is not a complete graph because the isolated vertex is not connected to the others.")
    print(f"The chromatic number of K_{n-1} is n - 1 = {n-1}.")
    print("The chromatic number of an isolated vertex is 1.")
    print(f"The chromatic number of G' is the maximum of its components, which is max({n-1}, 1) = {n-1}.\n")
    
    print("Step 4: Conclusion.")
    print(f"We have shown that any such graph G requires at most {n-1} colors, and that there exists at least one such graph that requires exactly {n-1} colors.")
    print("Therefore, the maximum number of colors we need is n - 1.\n")
    
    print("Final Calculation:")
    max_colors = n - 1
    print(f"Maximum colors = {n} - 1 = {max_colors}")

solve_coloring_problem()