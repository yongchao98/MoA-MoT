def solve_polyhedron_vertices():
    """
    Solves the mathematical problem about the possible number of vertices of a convex polyhedron
    with three quadrilateral projections onto planes in general position.

    The problem is a geometric one and the solution comes from mathematical reasoning,
    not computation. The possible numbers of vertices for such polyhedra are known to be 6, 8, and 12.

    - V=6: An octahedron. Projections along its 3 axes of 2-fold symmetry are rhombi (quadrilaterals).
      These axes are orthogonal and thus define planes in general position.
    - V=8: A cuboid (or any right quadrilateral prism). Projections along its 3 principal axes are
      rectangles (quadrilaterals). These axes are orthogonal.
    - V=12: A cuboctahedron. Projections along its 3 axes of 4-fold symmetry (passing through centers
      of opposite square faces) are squares (quadrilaterals). These axes are orthogonal.

    The code below simply prints these determined values.
    """
    possible_vertices = [6, 8, 12]
    print("The set of possible numbers of vertices is:")
    for v in possible_vertices:
        print(v)

solve_polyhedron_vertices()