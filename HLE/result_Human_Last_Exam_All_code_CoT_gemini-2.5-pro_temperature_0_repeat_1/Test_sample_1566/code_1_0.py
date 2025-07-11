def solve_topology_problem():
    """
    This function provides a step-by-step deduction to solve the topology problem
    and prints the final answer.
    """

    # The problem asks for the number of topologically distinct continua X with two properties:
    # 1. X has n end points, where 1 < n < infinity.
    # 2. X has exactly two orbits under its auto-homeomorphism group.

    # Step 1: Analyze the properties and their implications.
    # The two orbits must correspond to the set of end points (E) and the set of non-end-points (I).
    # This means all non-end-points in I are topologically equivalent to each other.

    # Step 2: Deduce the structure of the continuum.
    # The topological equivalence of all non-end-points is a very strong condition.
    # In a graph-like continuum, it means there can be no "junctions" or "branch points"
    # (vertices of degree > 2), as these would be topologically distinct from points
    # on an edge (which have degree 2).
    # This restricts the continuum to have a maximum degree of 2 everywhere.

    # Step 3: Identify possible structures.
    # A connected continuum with maximum degree 2 must be either:
    # a) A simple closed curve (homeomorphic to a circle).
    # b) A simple arc (homeomorphic to a line segment [0, 1]).

    # Step 4: Check the candidates against the given properties.
    
    # Candidate A: Simple Closed Curve (Circle)
    # - Number of end points: 0. This fails Property (1) which requires more than one end point.
    is_circle_a_solution = False

    # Candidate B: Simple Arc ([0, 1])
    # - Number of end points: 2. This satisfies Property (1) since 2 > 1.
    # - Number of orbits: 2. The two end points {0, 1} form one orbit, and the
    #   interior (0, 1) forms the second orbit. This satisfies Property (2).
    is_arc_a_solution = True

    # Step 5: Count the number of distinct solutions.
    # The only structure that works is the simple arc. Any continuum satisfying the
    # conditions must be homeomorphic to a simple arc. Therefore, all such continua
    # belong to the same topological equivalence class.
    
    number_of_solutions = 1 if is_arc_a_solution and not is_circle_a_solution else 0

    print("Based on the analysis, the only continuum that satisfies the given properties is a simple arc.")
    print("All simple arcs are topologically equivalent (homeomorphic).")
    print(f"Therefore, the number of topologically distinct continua is:")
    
    # Final equation output
    print(f"{number_of_solutions}")

solve_topology_problem()