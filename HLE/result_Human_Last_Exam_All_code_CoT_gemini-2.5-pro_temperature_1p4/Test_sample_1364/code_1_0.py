def solve_polyhedron_vertices():
    """
    This function determines the set of possible numbers of vertices for a convex polyhedron
    that has quadrilateral projections on three planes in general position.
    
    Based on geometric analysis, the polyhedra that satisfy these conditions are:
    1. The Tetrahedron (V=4): Projections along directions connecting the midpoints
       of opposite edges are quadrilaterals. The three such directions are mutually
       orthogonal for a regular tetrahedron.
    2. The Octahedron (V=6): Projections along its three main axes are squares.
       These axes are mutually orthogonal.
    3. The Cube (V=8): Projections along its face diagonals (e.g., directions
       (1,1,0), (1,-1,0), (0,1,1)) are rectangles. These directions can be chosen
       to be in a general position (linearly independent).
       
    Other polyhedra, like prisms or pyramids, fail because the directions that yield
    quadrilateral projections are not in general position (e.g., they are coplanar).
    """
    
    # The set of possible numbers of vertices.
    possible_vertices = {4, 6, 8}
    
    print("The set of possible numbers of vertices for such a polyhedron P is:")
    # We are asked to output each number in the final equation.
    # We will print the elements of the set.
    print(sorted(list(possible_vertices)))

solve_polyhedron_vertices()