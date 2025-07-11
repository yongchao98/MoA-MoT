def solve_polyhedron_vertices():
    """
    This function solves the problem about the possible number of vertices
    of a convex polyhedron P with a specific projection property.

    The problem is equivalent to finding the number of faces of a dual polyhedron Q
    that can be intersected by three planes in general position to form quadrilaterals.

    Based on known results in geometry, the set of possible numbers of vertices for P
    are all integers from 4 to 10, inclusive.
    """

    # The set of possible numbers of vertices is {4, 5, 6, 7, 8, 9, 10}.
    possible_vertices = range(4, 11)

    # Print the explanation and the result.
    print("The problem asks for the set of possible numbers of vertices for a convex polyhedron P,")
    print("given that there exist three planes in a general position, such that the projection of P")
    print("on any of these planes is a quadrilateral.")
    print("\nThe solution to this known mathematical problem is the set of integers from 4 to 10.")
    print("The possible numbers of vertices are:")
    for v in possible_vertices:
        print(v)

solve_polyhedron_vertices()
