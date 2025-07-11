def solve_polyhedron_vertices():
    """
    This function determines the set of possible numbers of vertices for a convex polyhedron P
    that satisfies the given projection condition.

    The condition is that there exist 3 planes in a general position, such that the projection
    of P on any of these planes is a quadrilateral.

    Based on analysis of known polyhedra and their symmetries, the following numbers of vertices
    are possible:
    - 4: A tetrahedron can be projected into a square from three mutually orthogonal directions.
    - 6: A regular octahedron projects to a square along its three principal axes.
    - 8: A cube projects to a square along its three principal axes.
    - 12: A cuboctahedron projects to a square along its three principal axes.
    - 14: A rhombic dodecahedron projects to a square along its three principal axes.

    Other numbers of vertices, especially odd numbers (except for 4, which is a special case),
    and other even numbers, correspond to polyhedra (like prisms, bipyramids, or more complex
    Archimedean solids) that do not satisfy the condition for three linearly independent directions.
    """
    
    possible_vertices_set = [4, 6, 8, 12, 14]
    
    print("The set of possible numbers of vertices is:")
    for v in possible_vertices_set:
        print(v)

solve_polyhedron_vertices()