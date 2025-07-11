def solve_robot_arm_problem():
    """
    This function determines the dimensions of the disjoint connected manifold
    components of the configuration space X_4 for a 4-segment robot arm.

    The derivation proceeds by stratifying the space X_4 based on the dimension
    of the vector space spanned by the four unit vectors (v1, v2, v3, v4)
    that define a configuration. The sum of these vectors must be zero.
    """

    # Dimension of the generic, non-planar component (span(v_i) is 3D).
    # The initial space is (S^2)^4, with dimension 4 * 2 = 8.
    # The constraint v1+v2+v3+v4=0 removes 3 degrees of freedom.
    # Resulting dimension: 8 - 3 = 5.
    # This forms a single connected component.
    dim_generic = 5
    y_1 = dim_generic
    print(f"The dimension of the largest component (non-planar configurations) is y_1 = {y_1}")

    # Dimension of the planar, non-collinear components (span(v_i) is 2D).
    # These configurations must form rhombi, implying pairwise opposite vectors.
    # There are 3 ways to pair up the 4 vectors.
    # Dimension: 2 (for the plane) + 1 (for v1) + 1 (for v2) = 4.
    # This results in 3 disjoint connected components of dimension 4.
    dim_planar = 4
    y_2, y_3, y_4 = dim_planar, dim_planar, dim_planar
    print(f"The dimensions of the three planar components are y_2 = {y_2}, y_3 = {y_3}, y_4 = {y_4}")


    # Dimension of the collinear components (span(v_i) is 1D).
    # This requires two vectors to be 'u' and two to be '-u'.
    # There are 3 ways to partition the vectors this way.
    # Each component is defined by the choice of 'u' on S^2, so the dimension is 2.
    # This results in 3 disjoint connected components of dimension 2.
    dim_collinear = 2
    y_5, y_6, y_7 = dim_collinear, dim_collinear, dim_collinear
    print(f"The dimensions of the three collinear components are y_5 = {y_5}, y_6 = {y_6}, y_7 = {y_7}")


    # The full list of dimensions, sorted in descending order.
    dimensions = (y_1, y_2, y_3, y_4, y_5, y_6, y_7)
    
    # Final formatted output
    final_answer_string = ",".join(map(str, dimensions))
    print("\nThe list of dimensions (y_1, ..., y_ell) is:")
    print(final_answer_string)

solve_robot_arm_problem()