def solve_polyhedron_vertices():
    """
    This function explains the reasoning and prints the set of possible numbers of vertices
    for a 3D convex polyhedron P with the given projection property.
    """

    print("Analyzing the possible number of vertices (V) for the polyhedron P.")
    print("-" * 60)

    # Case V=4
    print("Case V = 4 (e.g., a tetrahedron):")
    print("The projection of a tetrahedron from any generic viewpoint is a quadrilateral.")
    print("Therefore, we can easily choose 3 linearly independent generic directions.")
    print("Result: V = 4 is a possible number.\n")

    # Case V=5
    print("Case V = 5 (e.g., a square pyramid or triangular bipyramid):")
    print("A square pyramid has only one set of directions that give a quadrilateral projection (those nearly normal to the base). These directions are not linearly independent.")
    print("A triangular bipyramid can be projected into a quadrilateral from directions in its 'equatorial plane', but these directions are coplanar.")
    print("Projections from other directions typically yield triangles or pentagons.")
    print("Result: V = 5 is not a possible number.\n")

    # Case V is even and V >= 6
    print("Case V is an even number and V >= 6:")
    print("These can be constructed using n-gonal prisms, where V = 2n for n >= 3.")
    print("A projection nearly normal to any of the n quadrilateral side-faces is a quadrilateral.")
    print("We can choose three such faces whose normals are not collinear (for n>=3) and slightly perturb these normal vectors to make them linearly independent, while preserving the quadrilateral projection.")
    print("Result: All even numbers V >= 6 are possible (V = 6, 8, 10, ...).\n")

    # Case V is odd and V >= 7
    print("Case V is an odd number and V >= 7:")
    print("These can be constructed by 'capping' an n-gonal prism (V=2n) on one of its bases, creating a polyhedron with V = 2n + 1 vertices (for n>=3).")
    print("The original quadrilateral side-faces of the prism remain. The argument for even V still holds, as the new 'cap' does not obstruct the view of these faces.")
    print("Result: All odd numbers V >= 7 are possible (V = 7, 9, 11, ...).\n")
    print("-" * 60)

    # Final Answer
    print("CONCLUSION:")
    print("The set of possible numbers of vertices for such a polyhedron P is:")
    print("The number 4, and all integers greater than or equal to 6.")
    final_set_representation = "{4} U {V in Z | V >= 6}"
    print(f"In set notation: {final_set_representation}")


solve_polyhedron_vertices()