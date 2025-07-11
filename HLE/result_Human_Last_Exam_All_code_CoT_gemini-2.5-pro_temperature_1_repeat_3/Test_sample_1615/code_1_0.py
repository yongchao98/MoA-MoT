# The user wants to find the maximum number of colors needed to properly color a graph G with n=12345 vertices, given that G is not a complete graph.

# Plan:
# 1. Define the number of vertices, n.
# 2. Explain the fundamental upper bound for the chromatic number of a graph.
# 3. Explain how this bound changes with the constraint that the graph is not complete.
# 4. Demonstrate that this new, lower bound is achievable by constructing a valid example graph.
# 5. Conclude with the maximum number of colors and perform the calculation.

def solve_graph_coloring_problem():
    """
    Calculates and explains the maximum number of colors for a non-complete graph.
    """
    # 1. Define the number of vertices.
    num_vertices = 12345

    # 2. Explain the general principle.
    print("Let n be the number of vertices in a graph G.")
    print(f"In this problem, n = {num_vertices}.")
    print("\nA fundamental result in graph theory states that the chromatic number of any graph G, denoted χ(G), is less than or equal to its number of vertices n.")
    print("χ(G) ≤ n")
    print("\nThe only graph that requires exactly n colors is the complete graph K_n, where every vertex is connected to every other vertex.")

    # 3. Apply the problem's specific constraint.
    print("\nThe problem specifies that our graph G is NOT a complete graph.")
    print("This means there must be at least one pair of vertices, say u and v, that are not adjacent (not connected by an edge).")
    print("Because u and v are not adjacent, they can be assigned the same color in a proper coloring.")
    print("This fact immediately establishes a tighter upper bound for any non-complete graph: χ(G) ≤ n - 1.")

    # 4. Show that the upper bound is achievable.
    print("\nTo determine if n - 1 is the maximum possible number of colors, we must check if this bound is achievable.")
    print("Consider a graph constructed by taking a complete subgraph on n-1 vertices (K_{n-1}) and adding one additional, isolated vertex.")
    print("This graph has n vertices but is not a complete graph.")
    print(f"The K_{{n-1}} subgraph by itself requires n-1 distinct colors.")
    print("Since this graph exists and requires n-1 colors, the maximum is indeed n-1.")

    # 5. Perform the final calculation.
    max_colors = num_vertices - 1
    print("\nTherefore, the maximum number of colors needed is n - 1.")
    print("\nFinal Calculation:")
    print(f"{num_vertices} - 1 = {max_colors}")

solve_graph_coloring_problem()