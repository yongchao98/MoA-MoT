def solve_great_icosahedron_triangulation():
    """
    Calculates the number of triangular faces on the convex hull of a great icosahedron.

    The convex hull of a great icosahedron is a regular icosahedron. This function
    calculates the number of faces of a regular icosahedron from its vertex and edge properties.
    """
    # A regular icosahedron has 12 vertices.
    vertices = 12

    # At each vertex of a regular icosahedron, 5 edges meet.
    edges_per_vertex = 5

    # To calculate the total number of edges (E), we multiply vertices by edges_per_vertex
    # and divide by 2 (since each edge connects two vertices).
    # Equation: E = (V * 5) / 2
    total_edges = (vertices * edges_per_vertex) / 2

    # The faces of an icosahedron are triangles, so each face has 3 edges.
    edges_per_face = 3

    # To calculate the total number of faces (F), we multiply the total edges by 2
    # (since each edge is shared by two faces) and divide by the edges per face.
    # Equation: F = (E * 2) / 3
    total_faces = (total_edges * 2) / edges_per_face

    print("The outer-hull is a regular icosahedron.")
    print(f"Step 1: Calculate the number of edges from {vertices} vertices, with {edges_per_vertex} edges per vertex.")
    print(f"Equation: ({vertices} * {edges_per_vertex}) / 2 = {int(total_edges)} edges")
    print("")
    print(f"Step 2: Calculate the number of faces from {int(total_edges)} edges, with each face being a triangle ({edges_per_face} edges).")
    print(f"Equation: ({int(total_edges)} * 2) / {edges_per_face} = {int(total_faces)} triangles")
    print("")
    print(f"The minimal number of triangles is {int(total_faces)}.")

solve_great_icosahedron_triangulation()
<<<20>>>