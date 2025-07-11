def solve_great_icosahedron_triangles():
    """
    Calculates the number of visible triangles on the surface of a great icosahedron.
    This corresponds to the number of faces of its convex hull, which is a regular icosahedron.
    """
    # A regular icosahedron has 12 vertices and 30 edges.
    vertices = 12
    edges = 30

    # According to Euler's formula for polyhedra: V - E + F = 2.
    # We can find the number of faces (F) by rearranging the formula to F = 2 - V + E.
    # Since the faces of an icosahedron are triangles, this number is our answer.
    num_faces = 2 - vertices + edges

    # The problem asks to output the numbers in the final equation.
    print(f"The outer-hull of the great icosahedron is a regular icosahedron.")
    print(f"The number of vertices (V) is {vertices}.")
    print(f"The number of edges (E) is {edges}.")
    print(f"Using Euler's formula V - E + F = 2, we calculate the number of faces (F).")
    print(f"The final equation is: F = 2 - {vertices} + {edges}")
    print(f"The minimal number of visible triangles is: {num_faces}")

solve_great_icosahedron_triangles()