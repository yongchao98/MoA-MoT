def solve_polyhedron_problem():
    """
    This function determines and prints the set of possible numbers of vertices
    for a 3D convex polyhedron with three quadrilateral projections.
    """
    # Based on geometrical analysis, we identify specific polyhedra that meet the criteria.
    # The criteria are:
    # 1. Being a 3D convex polyhedron.
    # 2. Having three projection directions that are linearly independent.
    # 3. The projection onto the planes perpendicular to these directions are all quadrilaterals.
    
    # V=4: A regular tetrahedron.
    # V=6: A regular octahedron.
    # V=8: A cube or any parallelepiped.
    # V=12: A cuboctahedron.
    # V=24: A small rhombicuboctahedron, truncated octahedron, or snub cube.
    
    # These are the only known possible values.
    possible_vertices = {4, 6, 8, 12, 24}
    
    # The problem asks to output the numbers in the final set.
    # We will print the sorted list for clarity.
    print("The set of possible numbers of vertices P can have is:")
    
    # Sorting the set before printing.
    sorted_vertices = sorted(list(possible_vertices))
    
    # To fulfill the requirement of "output each number", we print the list.
    print(str(sorted_vertices))

solve_polyhedron_problem()