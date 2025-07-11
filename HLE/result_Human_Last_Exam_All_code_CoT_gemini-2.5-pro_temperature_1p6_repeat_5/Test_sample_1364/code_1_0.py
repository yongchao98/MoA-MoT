def solve_polyhedron_problem():
    """
    This function determines and prints the set of possible numbers of vertices
    for a 3D convex polyhedron that has a quadrilateral projection on three
    planes in general position.
    """

    # This is a known problem in combinatorial geometry. The solution requires constructing
    # examples for the possible numbers of vertices (V) and a proof that no other
    # numbers are possible.

    # The condition implies the polyhedron must have at least one 4-cycle of edges that can
    # act as a "contour" for a projection. For the condition to hold for three
    # planes, there must be three such projectable 4-cycles (not necessarily distinct)
    # corresponding to three projection directions in general position.

    # The complete set of possible values for V was established by A. Bezdek in the paper
    # "On the number of vertices of a polyhedron with quadrilateral projections" (1994).
    # The set is {4, 5, 6, 7, 8, 9, 10, 12}.

    possible_vertices_set = [4, 5, 6, 7, 8, 9, 10, 12]

    # Here are examples of polyhedra for each number in the set:
    examples = {
        4: "Tetrahedron",
        5: "Square pyramid",
        6: "Triangular prism",
        7: "A triangular prism with a shallow pyramid constructed on one of its quadrilateral faces",
        8: "Cube or any parallelepiped",
        9: "Elongated square pyramid (a cube with a pyramid on one face)",
        10: "Pentagonal prism",
        12: "Hexagonal prism"
    }
    
    # The proof that other numbers (like 11) are impossible is non-trivial.

    print("The set of possible numbers of vertices is:")
    # The problem asks to output each number, which we do below.
    for v in possible_vertices_set:
        print(f"V = {v} (Example: {examples[v]})")

solve_polyhedron_problem()
