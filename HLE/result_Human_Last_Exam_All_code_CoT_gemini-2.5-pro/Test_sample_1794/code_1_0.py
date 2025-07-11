def solve_feynman_diagram_problem():
    """
    This script solves a two-part problem about 3-loop Feynman diagrams
    in scalar phi^3 theory.
    """

    # --- Part 1: Counting Planar Graphs ---
    print("--- Solving Part 1 ---")
    print("Counting the number of planar 3-loop 4-point graphs without 'vertex corrections'.")
    print("This is interpreted as counting the 1PI graphs that do not contain a triangle subgraph.")

    # The 6 known 1PI planar diagrams can be classified by this property.
    graphs_with_triangles = [
        "Twisted Tennis Court", 
        "Penta-box variant 1", 
        "Penta-box variant 2"
    ]
    graphs_without_triangles = [
        "Ladder", 
        "Tennis Court", 
        "Zig-zag"
    ]
    
    num_without_triangles = len(graphs_without_triangles)
    
    print("\nThe diagrams without triangle subgraphs are:")
    for graph_name in graphs_without_triangles:
        print(f"- {graph_name}")

    # Displaying the calculation as a sum
    equation_str = " + ".join(["1"] * num_without_triangles)
    
    print(f"\nThe total count is the sum of these graphs: {equation_str} = {num_without_triangles}")
    print(f"Result for Part 1: The number of distinct planar graphs is {num_without_triangles}.\n")


    # --- Part 2: Power of the Leading Divergence ---
    print("--- Solving Part 2 ---")
    print("Finding the power of the leading divergent term of the Feynman integral's epsilon expansion.")
    
    # For massless on-shell diagrams, divergences are of Infrared (IR) origin.
    # The leading pole 1/epsilon^k in an L-loop planar amplitude has a maximal power k = 2*L.
    
    # Number of loops for the diagram
    L = 3
    
    # Calculate the power of the leading pole
    k = 2 * L
    
    print(f"The calculation is for L = {L}-loop diagrams.")
    print("The power 'k' of the leading divergent term (1/epsilon^k) is given by the formula: k = 2 * L.")
    print(f"Plugging in L={L}, we get k = 2 * {L} = {k}")
    
    print(f"Result for Part 2: The number of the power of the leading divergent term is {k}.")


# Run the solver function
solve_feynman_diagram_problem()