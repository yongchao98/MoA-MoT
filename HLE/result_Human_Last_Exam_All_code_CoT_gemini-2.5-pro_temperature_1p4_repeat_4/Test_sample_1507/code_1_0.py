def solve_triangulation():
    """
    Calculates the minimal number of triangles needed for triangulating
    the outer-hull of the great icosahedron.
    """
    
    # Step 1: The outer hull of a great icosahedron is a regular icosahedron.
    # A regular icosahedron is a convex polyhedron.
    
    # Step 2: The surface of a regular icosahedron is already composed of triangles.
    # The number of triangles in a minimal triangulation is therefore its number of faces.
    
    # Step 3: A regular icosahedron has 20 faces.
    num_faces_of_icosahedron = 20
    
    print("The outer hull of a great icosahedron is a regular icosahedron.")
    print("The surface of a regular icosahedron is already triangulated by its own faces.")
    print("Therefore, the minimal number of triangles is the number of faces of the icosahedron.")
    
    # Final equation and result
    print("\nFinal Equation:")
    print(f"Number of Triangles = Number of Faces")
    print(f"Number of Triangles = {num_faces_of_icosahedron}")

solve_triangulation()