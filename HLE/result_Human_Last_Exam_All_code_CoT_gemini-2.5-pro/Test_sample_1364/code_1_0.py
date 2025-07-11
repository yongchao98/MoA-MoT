def solve():
    """
    This function determines the set of possible numbers of vertices for a convex polyhedron P,
    given that P can be projected onto three planes in general position as a quadrilateral.

    Based on geometric analysis, we test several classes of polyhedra:
    1.  Tetrahedron (4 vertices): A regular tetrahedron can be projected as a square onto the
        three coordinate planes. A square is a quadrilateral. So, 4 is a possible number.
    2.  Octahedron (6 vertices): A regular octahedron can be projected as a square onto the
        three coordinate planes. So, 6 is a possible number.
    3.  Cube (8 vertices): A cube can be projected as a square onto the three coordinate
        planes. So, 8 is a possible number.
    4.  Other polyhedra like prisms (except the cube), pyramids, and bipyramids can be shown
        to fail the condition. For a pyramid with a quadrilateral base, only one projection
        (onto the base) is a quadrilateral; others are triangles. For a prism, all directions
        yielding a quadrilateral projection are coplanar, failing the "general position" criteria.

    Therefore, the set of possible numbers of vertices is {4, 6, 8}.
    """
    possible_vertices = {4, 6, 8}
    
    # The problem requires printing the final equation, showing each number.
    # We will print the elements of the set.
    print("The set of possible numbers of vertices is {", end="")
    
    # Sort the numbers for a consistent output format
    sorted_vertices = sorted(list(possible_vertices))
    
    for i, v in enumerate(sorted_vertices):
        if i > 0:
            print(", ", end="")
        print(f"{v}", end="")
        
    print("}")

solve()