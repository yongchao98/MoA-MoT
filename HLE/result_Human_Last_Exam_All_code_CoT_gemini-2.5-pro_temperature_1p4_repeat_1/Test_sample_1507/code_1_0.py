def solve_triangulation():
    """
    Calculates the number of triangles for triangulating the outer-hull
    of the great icosahedron.
    """
    # The outer hull of a great icosahedron is a regular icosahedron.
    # We can use Euler's formula (V - E + F = 2) for convex polyhedra to find the number of faces (F).
    # For a regular icosahedron:
    vertices = 12
    edges = 30

    # The equation is: V - E + F = 2
    # We solve for F: F = 2 - V + E
    faces = 2 - vertices + edges

    print("The outer hull of the great icosahedron is a regular icosahedron.")
    print("The number of visible triangles equals the number of faces of this icosahedron.")
    print("We use Euler's formula: V - E + F = 2")
    print(f"A regular icosahedron has {vertices} vertices (V) and {edges} edges (E).")
    print("\nCalculating the number of faces (F):")
    print("-----------------------------------")
    print(f"Equation: {vertices} - {edges} + F = 2")
    print(f"Solving for F: F = 2 - {vertices} + {edges}")
    print(f"Final calculation: F = {faces}")
    print("-----------------------------------")
    print(f"\nThe minimal number of triangles needed is {faces}.")

solve_triangulation()