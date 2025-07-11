def solve():
    """
    This function determines the set of possible numbers of vertices for a convex polyhedron
    that can be projected onto three planes in general position to form quadrilaterals.
    
    The reasoning is based on constructing such polyhedra:
    1.  V=6: A regular octahedron works. Projections along its axes of symmetry are squares.
    2.  V=8+2k: A cube (k=0) and its truncations (k=1 to 8) work. Projections along
        directions like (1,1,0), (1,0,1), and (0,1,1) are rectangles.
    
    This leads to the following set of possible vertex counts.
    """
    
    # Construction for V=6 (Octahedron)
    v_octahedron = 6
    
    # Construction for V=8+2k (Cube and its truncations)
    v_cube_family = []
    for k in range(0, 9): # k from 0 to 8 truncations
        v_cube_family.append(8 + 2 * k)
        
    possible_vertices = sorted(list(set([v_octahedron] + v_cube_family)))
    
    print("Based on geometric constructions, the set of possible numbers of vertices is:")
    print(possible_vertices)

solve()