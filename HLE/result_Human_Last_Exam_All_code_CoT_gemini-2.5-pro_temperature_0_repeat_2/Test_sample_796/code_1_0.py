def solve_robot_arm_problem():
    """
    This function determines the dimensions of the manifolds that form
    a minimal disjoint decomposition of the configuration space X_4.

    The configuration space X_4 of 4-segment unit-length robot arms
    that form a closed loop is stratified by the length of the diagonal
    vector connecting the start of the first segment to the end of the second.
    This length, r = ||v1 + v2||, can range from 0 to 2.

    This stratification yields a disjoint union of three connected manifolds:
    1.  r in (0, 2): The generic case. This is a 5-dimensional manifold.
        dim = dim((S^2)^4) - dim(R^3) = 4*2 - 3 = 5.
    2.  r = 0: This implies v2 = -v1 and v4 = -v3. The space is parameterized
        by (v1, v3) in S^2 x S^2. This is a 4-dimensional manifold.
        dim = dim(S^2) + dim(S^2) = 2 + 2 = 4.
    3.  r = 2: This implies v2 = v1 and v4 = v3. The constraint v1+v2+v3+v4=0
        becomes 2*v1 + 2*v3 = 0, so v3 = -v1. The space is parameterized
        by v1 in S^2. This is a 2-dimensional manifold.
        dim = dim(S^2) = 2.

    The dimensions, sorted in descending order, are (5, 4, 2).
    """
    y1 = 5
    y2 = 4
    y3 = 2
    
    # The problem asks to output each number in the final equation.
    # We interpret this as printing the components of the dimension tuple.
    print(f"{y1},{y2},{y3}")

solve_robot_arm_problem()