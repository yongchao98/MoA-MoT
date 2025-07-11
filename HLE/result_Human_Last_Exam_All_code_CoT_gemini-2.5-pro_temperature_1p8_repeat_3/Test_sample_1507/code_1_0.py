def solve_great_icosahedron_triangulation():
    """
    Calculates the number of triangles on the outer-hull of a great icosahedron.

    The outer-hull of a great icosahedron is a regular icosahedron.
    The number of triangles needed to triangulate this hull is simply its number of faces.
    We use Euler's formula (V - E + F = 2) to find the number of faces (F).
    """

    # Properties of a regular icosahedron (the outer-hull)
    V = 12  # Number of vertices
    E = 30  # Number of edges
    Euler_char = 2 # Euler characteristic for convex polyhedra

    # Using Euler's formula to find the number of faces (F): F = 2 - V + E
    F = Euler_char - V + E

    print("The outer-hull of the great icosahedron is a regular icosahedron.")
    print("The number of visible triangles is equal to the number of faces (F) of this icosahedron.")
    print("We can find F using Euler's formula: V - E + F = 2, or F = 2 - V + E.")
    print(f"Given Vertices (V) = {V} and Edges (E) = {E}:")
    # As requested, outputting each number in the final equation
    print(f"The calculation is: F = {Euler_char} - {V} + {E}")
    print(f"So, the number of triangles is {F}.")

solve_great_icosahedron_triangulation()
<<<20>>>