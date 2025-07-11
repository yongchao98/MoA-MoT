def calculate_icosahedron_faces():
    """
    Calculates the number of faces of a regular icosahedron, which corresponds to
    the number of triangles on the outer-hull of a great icosahedron.
    """
    # The outer-hull (convex hull) of a great icosahedron is a regular icosahedron.
    # We need to find the number of faces of a regular icosahedron.

    # Properties of a regular icosahedron:
    num_vertices = 12  # V
    num_edges = 30     # E

    # Euler's formula for polyhedra is V - E + F = 2, where F is the number of faces.
    # We can rearrange it to find F: F = 2 - V + E.
    num_faces = 2 - num_vertices + num_edges

    print("The outer-hull of the great icosahedron is a regular icosahedron.")
    print("A regular icosahedron has {} vertices (V) and {} edges (E).".format(num_vertices, num_edges))
    print("\nUsing Euler's formula for polyhedra (V - E + F = 2), we can find the number of faces (F).")
    print("The equation is: F = 2 - V + E")
    print("Substituting the values:")
    # The prompt requires outputting each number in the final equation.
    print("F = {} - {} + {}".format(2, num_vertices, num_edges))
    print("F = {}".format(num_faces))
    
    print("\nThe faces of a regular icosahedron are triangles.")
    print("Therefore, the minimal number of triangles needed for triangulating its surface is equal to its number of faces.")
    print("\nMinimal number of triangles = {}".format(num_faces))

calculate_icosahedron_faces()