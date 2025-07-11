def solve_great_icosahedron_triangles():
    """
    Calculates the number of visible triangles on the outer-hull of a great icosahedron.

    The outer-hull (or convex hull) of the great icosahedron is a regular icosahedron.
    A regular icosahedron is a Platonic solid that has 20 faces.
    Each face of a regular icosahedron is an equilateral triangle.

    Therefore, the number of visible triangles from the outside is equal to the
    number of faces of a regular icosahedron.
    """
    
    # Number of faces of a regular icosahedron
    num_faces = 20
    
    print("The outer-hull of the great icosahedron is a regular icosahedron.")
    print("The number of faces (which are triangles) on a regular icosahedron is:")
    print(num_faces)

solve_great_icosahedron_triangles()