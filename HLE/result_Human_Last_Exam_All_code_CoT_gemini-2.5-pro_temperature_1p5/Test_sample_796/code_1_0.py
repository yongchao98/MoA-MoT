def solve_robot_arm_dimensions():
    """
    This function determines the dimensions of the connected manifold components
    of the configuration space X_4.

    The space X_4 is stratified based on the geometric properties of the vector arrangements.
    The strata are:
    1. The set of non-coplanar configurations, which is a manifold of dimension 5.
    2. The set of coplanar but not collinear configurations, a manifold of dimension 4.
    3. The set of collinear configurations, which splits into 3 disjoint manifolds,
       each of dimension 2.

    The dimensions of these components, sorted in descending order, are calculated.
    """

    # The dimensions of the disjoint connected manifold components
    y1 = 5  # Dimension of the manifold of non-coplanar configurations
    y2 = 4  # Dimension of the manifold of coplanar, non-collinear configurations
    y3 = 2  # Dimension of the first type of collinear configurations
    y4 = 2  # Dimension of the second type of collinear configurations
    y5 = 2  # Dimension of the third type of collinear configurations

    # The list of dimensions in descending order
    dimensions = [y1, y2, y3, y4, y5]

    # Print the result in the specified format
    print(",".join(map(str, dimensions)))

solve_robot_arm_dimensions()