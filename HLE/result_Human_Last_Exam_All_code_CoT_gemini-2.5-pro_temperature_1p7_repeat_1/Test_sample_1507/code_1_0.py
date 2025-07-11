def solve_triangulation():
    """
    Calculates the minimal number of triangles for triangulating the
    outer-hull of the great icosahedron.
    """
    # The outer-hull of a great icosahedron is a regular icosahedron.
    # A regular icosahedron is a Platonic solid with 20 faces.
    # All of its faces are equilateral triangles.
    num_faces_regular_icosahedron = 20

    # Since the surface of a regular icosahedron is already composed of triangles,
    # the minimal number of triangles for its triangulation is simply its number of faces.
    print("The outer-hull of the great icosahedron is a regular icosahedron.")
    print("A regular icosahedron is already composed of triangular faces.")
    print("Therefore, the minimal number of triangles needed is equal to its number of faces.")
    print("\nFinal Calculation:")
    # The prompt requires printing the numbers in the final equation.
    # We will format the output to look like a statement or equation.
    equation_part_1 = "Minimal number of triangles"
    equation_part_2 = str(num_faces_regular_icosahedron)
    print(f"{equation_part_1} = {equation_part_2}")

solve_triangulation()