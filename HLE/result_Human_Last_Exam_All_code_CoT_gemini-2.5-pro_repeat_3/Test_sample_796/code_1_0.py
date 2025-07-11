def solve_robot_arm_problem():
    """
    This function determines the dimensions of the manifold components of the robot arm configuration space.

    The configuration space X_4 of a 4-segment unit-length robot arm is the set of vectors
    (v1, v2, v3, v4) in (S^2)^4 such that v1 + v2 + v3 + v4 = 0.

    This space can be decomposed into a disjoint union of connected manifolds without boundary.
    This decomposition is given by the strata of the space, defined by the rank of the
    derivative of the constraint function.

    1. Generic Stratum (rank 3): The set of non-collinear configurations.
       The dimension is dim((S^2)^4) - dim(R^3) = (4 * 2) - 3 = 5.
       This forms one connected component.

    2. Singular Stratum (rank 2): The set of collinear configurations.
       For the vectors to be collinear and sum to zero, two must be a unit vector `u`
       and two must be `-u`.
       The number of ways to assign these vectors to the four positions is C(4,2) = 6.
       Each of these 6 configuration types is parameterized by the choice of `u` from S^2,
       making each a manifold of dimension 2.

    Therefore, the space decomposes into one manifold of dimension 5 and six manifolds of dimension 2.
    The list of dimensions is (5, 2, 2, 2, 2, 2, 2).
    """
    
    # The dimensions of the disjoint manifold components
    y1 = 5
    y2 = 2
    y3 = 2
    y4 = 2
    y5 = 2
    y6 = 2
    y7 = 2
    
    # The problem asks for the dimensions in a specific format.
    # We print the result as a comma-separated string.
    print(f"{y1},{y2},{y3},{y4},{y5},{y6},{y7}")

solve_robot_arm_problem()