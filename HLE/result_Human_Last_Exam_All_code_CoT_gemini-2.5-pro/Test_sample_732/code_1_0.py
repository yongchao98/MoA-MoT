def solve_graph_problem():
    """
    Calculates and explains the maximum number of edges in a C4-free graph with 8 vertices.
    """
    n = 8
    forbidden_cycle_length = 4
    
    print(f"This problem asks for the maximum number of edges in a simple graph with {n} vertices that contains no quadrilaterals (C{forbidden_cycle_length}).")
    print("This value is known as the Turan number ex(8, C4).\n")
    
    # Step 1: Establish a lower bound
    print("--- Step 1: Establish a lower bound ---")
    print("We can construct a C4-free graph on 8 vertices to find a possible number of edges.")
    print("1. Start with the known maximal C4-free graph on 7 vertices, the Friendship Graph F_3.")
    print("   - F_3 consists of 3 triangles sharing a common central vertex.")
    print("   - It has 7 vertices and 9 edges.")
    print("2. Add an 8th vertex and connect it to the central vertex of F_3.")
    print("   - The new graph has 8 vertices.")
    print("   - The number of edges is 9 + 1 = 10.")
    print("   - This graph is C4-free because any new cycle would need to pass through the new vertex, which only has one neighbor.")
    lower_bound = 10
    print(f"This construction shows that the maximum number of edges is at least {lower_bound}.\n")
    
    # Step 2: Establish an upper bound
    print("--- Step 2: Establish an upper bound ---")
    print("We prove by contradiction that a C4-free graph with 8 vertices cannot have 11 edges.")
    print("1. Assume such a graph G exists, with n=8 vertices and m=11 edges.")
    print("2. The sum of degrees in G is 2 * m = 2 * 11 = 22.")
    print("3. The average degree is 22 / 8 = 2.75.")
    print("4. Therefore, there must be a vertex 'v' with degree at most 2.")
    print("5. Remove 'v'. The remaining graph G' has n'=7 vertices and m' = 11 - d(v) >= 11 - 2 = 9 edges.")
    print("6. We know ex(7, C4) = 9 (the friendship graph F_3). Since G' must be C4-free, m' must be exactly 9, and G' must be F_3.")
    print("7. This implies G is formed by adding a vertex of degree 2 to F_3.")
    print("8. However, adding a vertex 'v' connected to any two vertices 'x' and 'y' in F_3 will create a C4, because any two vertices in F_3 are either adjacent or have a common neighbor.")
    print("9. This is a contradiction. The assumption that a graph with 11 edges exists is false.\n")

    # Step 3: Conclusion
    print("--- Step 3: Conclusion ---")
    upper_bound = 11
    max_edges = 10
    print(f"We have shown that the maximum number of edges M satisfies M >= {lower_bound} and M < {upper_bound}.")
    print("Since the number of edges must be an integer, the maximum number of edges is 10.")
    
    final_equation_parts = {
        "vertices": n,
        "cycle_length": forbidden_cycle_length,
        "max_edges": max_edges
    }
    
    print("\nFinal equation: ex(n, Ck) = M")
    print(f"n = {final_equation_parts['vertices']}")
    print(f"k = {final_equation_parts['cycle_length']}")
    print(f"M = {final_equation_parts['max_edges']}")
    
solve_graph_problem()