def count_adjustable_graphs():
    """
    This function determines the number of non-isomorphic graphs
    with the specified properties by following a logical deduction.
    """
    
    # Properties of the graph G:
    # 1. Connected
    # 2. 3-regular
    # 3. 2000 vertices
    # 4. Has a perfect matching
    # 5. Is an "adjustable graph", meaning it has an adjustable perfect matching M.

    # An adjustable perfect matching M imposes a strong structure on G.
    # The non-matching edges of G form a 2-regular graph.
    # The connectivity of G implies this 2-regular graph must consist of two 1000-vertex components.
    # Let these components be C1 and C2.
    # The vertices of C1 must be the matching partners of the vertices of C2.
    # The adjustable property implies that if vertices v_i and v_j are adjacent in C1,
    # their partners u_i and u_j must be adjacent in C2.
    # This forces C1 and C2 to be isomorphic.
    
    # A connected 2-regular graph is a simple cycle.
    # Thus, C1 and C2 must both be 1000-cycles (C_1000).
    
    # The resulting graph G is constructed from two C_1000 cycles,
    # with a perfect matching connecting corresponding vertices.
    # This structure is the prism graph C_1000 x K_2.
    
    # All graphs with this structure are isomorphic to each other, regardless of labeling.
    # Therefore, there is only one such non-isomorphic graph.
    
    num_non_isomorphic_graphs = 1
    
    equation_explanation = "The number of non-isomorphic graphs is determined by the unique structure forced by the given properties."
    final_equation = f"Number of non-isomorphic graphs = {num_non_isomorphic_graphs}"
    
    print(equation_explanation)
    print(final_equation)

# Run the function to find the answer
count_adjustable_graphs()
<<<1>>>