def calculate_triangulation_faces():
    """
    Calculates the number of triangles in a minimal triangulation of the
    outer-hull of the great icosahedron.
    """
    # Step 1: The outer-hull of the great icosahedron is a regular icosahedron.
    # A regular icosahedron has 12 vertices.
    V = 12

    # Step 2: For any triangulation of a convex polyhedron, the number of faces (F)
    # can be determined from the number of vertices (V) using a formula derived
    # from Euler's characteristic (V - E + F = 2) and the edge-face relationship (3F = 2E).
    # The derived formula is: F = 2 * (V - 2).
    
    # Step 3: Calculate the number of triangles (F).
    F = 2 * (V - 2)
    
    print("The outer-hull of the great icosahedron is a regular icosahedron.")
    print(f"The number of vertices (V) for this shape is {V}.")
    print("\nThe number of triangles (F) can be calculated with the formula: F = 2 * (V - 2).")
    
    print("\nThe final equation with the numbers filled in is:")
    # We output the final equation showing all the numbers involved.
    print(f"{F} = 2 * ({V} - 2)")
    
    print(f"\nThus, the minimal number of triangles needed is {F}.")

calculate_triangulation_faces()