def solve_mobius_forest_problem():
    """
    Calculates the number of higher dimensional rooted forests on a standard
    triangulation of the Möbius band that fail to simplicially collapse onto their root.
    """

    # Step 1: Define the properties of the standard triangulation of the Möbius band.
    # This triangulation has f_0=6 vertices, f_1=12 edges, and f_2=6 triangles.
    # Its boundary is a cycle of 6 edges.
    num_triangles = 6
    total_edges = 12
    boundary_edges = 6
    interior_edges = total_edges - boundary_edges

    # Step 2: Apply the formula for 2-manifolds.
    # S(K) = 8 * (number of interior edges) + 4 * (number of triangles)
    term1_factor = 8
    term2_factor = 4
    
    result = term1_factor * interior_edges + term2_factor * num_triangles

    # Step 3: Print the final equation with the numbers.
    print(f"{term1_factor} * {interior_edges} + {term2_factor} * {num_triangles} = {result}")

solve_mobius_forest_problem()