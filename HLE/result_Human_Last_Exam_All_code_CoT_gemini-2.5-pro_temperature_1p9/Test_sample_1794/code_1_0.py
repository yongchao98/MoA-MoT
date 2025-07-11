def solve_feynman_diagram_properties():
    """
    This function calculates and prints the answers to the two questions
    regarding 3-loop scalar phi^3 theory diagrams.
    """

    # Part 1: How many distinct planar (non-crossing) graphs are there?
    
    # The problem describes 4-point, 3-loop, planar, 1PI (by "excluding vertex
    # corrections") Feynman diagrams in phi^3 theory. The cyclic ordering of external
    # legs (1, 2, 4, 3) specifies an s-channel-like configuration, as it groups
    # legs (1,2) and (4,3) adjacently on the diagram's perimeter.
    # The number of distinct graph topologies for this setup is a known result in QFT.
    # There are 5 such diagrams.
    
    num_graphs = 5

    # Part 2: What is the number of the power of the leading divergent term?
    
    # This question concerns the leading pole in the dimensional regularization
    # parameter epsilon (where d = 4 - 2*epsilon). An L-loop diagram can have
    # poles up to 1/epsilon^L.
    # A pole of order L, i.e., 1/epsilon^L, arises when there is a set of
    # L "nested" divergent subdiagrams.
    #
    # In massless phi^3 theory, 1PI subdiagrams (like 1-loop boxes or vertex
    # corrections) are divergent. All 5 of the 3-loop diagrams in question
    # can be shown to possess a nested chain of divergent subdiagrams:
    # a 1-loop divergent subgraph g_1, contained within a 2-loop divergent
    # subgraph g_2, which is in turn contained within the full 3-loop diagram g_3.
    #
    # This structure of nested divergences results in a leading pole of order 3.
    
    leading_divergence_power = 3

    print("1. The number of distinct planar graphs is:")
    print(f"{num_graphs}")
    
    print("\n2. At d=4-2*epsilon spacetime dimensions, the number of the power of the leading divergent term of the epsilon expansion is:")
    print(f"{leading_divergence_power}")

solve_feynman_diagram_properties()