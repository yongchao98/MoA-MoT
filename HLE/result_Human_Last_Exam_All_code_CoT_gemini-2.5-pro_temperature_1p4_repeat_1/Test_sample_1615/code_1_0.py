def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph
    with 12345 vertices that is not a complete graph.
    """
    
    # The number of vertices in the graph G
    n = 12345
    
    # The problem asks for the maximum possible chromatic number for a graph with n vertices
    # that is not a complete graph.
    
    print("Step 1: Establishing an upper bound for the chromatic number.")
    print(f"Let G be a graph with n = {n} vertices.")
    print("If G is not a complete graph, there must be at least two vertices, u and v, that are not adjacent.")
    print("We can assign the same color to u and v.")
    print(f"The remaining {n} - 2 vertices can each be assigned a unique new color.")
    print(f"This results in a proper coloring using at most 1 + ({n} - 2) = {n-1} colors.")
    print(f"Thus, the chromatic number chi(G) is at most {n-1}.")
    print("\n" + "="*50 + "\n")
    
    print("Step 2: Showing the upper bound is achievable.")
    print("Consider a graph G constructed from two disjoint components:")
    print(f"1. A complete graph K_{n-1}, which is a K_{{{n-1}}}.")
    print("2. A single isolated vertex K_1.")
    print(f"This graph G has ({n-1}) + 1 = {n} vertices and is not a complete graph.")
    print(f"The chromatic number of K_{{{n-1}}} is {n-1}.")
    print("The chromatic number of an isolated vertex is 1.")
    print(f"The chromatic number of G is the maximum of its components, which is {n-1}.")
    print("\n" + "="*50 + "\n")

    # Calculate the maximum number of colors
    max_colors = n - 1
    
    print("Conclusion:")
    print("The maximum number of colors needed for a proper vertex coloring of a graph with")
    print(f"{n} vertices, which is not a complete graph, is n - 1.")
    print("\nFinal Calculation:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()