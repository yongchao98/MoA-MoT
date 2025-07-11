def solve_critical_points_on_torus():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using its Betti numbers.
    """
    # The Betti numbers for the 2-torus (T^2) are well-known topological invariants.
    # b_0 is the number of connected components. The torus is one connected piece.
    b0 = 1
    
    # b_1 is the number of 1-dimensional "holes" or independent loops.
    # A torus has two such loops (one "around" the tube, one "through" the hole).
    b1 = 2
    
    # b_2 is the number of 2-dimensional "voids". The surface of the torus encloses one.
    b2 = 1
    
    # A result from Morse theory states that the number of critical points N of a
    # smooth function on a compact manifold is bounded below by the sum of its Betti numbers.
    # N >= b_0 + b_1 + b_2
    # This bound is sharp, meaning a function exists with exactly this many points.
    minimal_critical_points = b0 + b1 + b2
    
    # Print the equation as requested
    print(f"The minimal number of critical points N is bounded by the sum of Betti numbers:")
    print(f"N >= b_0 + b_1 + b_2")
    print(f"N >= {b0} + {b1} + {b2} = {minimal_critical_points}")
    print(f"\nThus, the minimal number of critical points is {minimal_critical_points}.")

solve_critical_points_on_torus()