def solve_configuration_space_dimensions():
    """
    This function determines the dimensions of the disjoint connected manifold
    components of the configuration space X_4.

    The space X_4 is stratified by the rank of the set of vectors {v1, v2, v3, v4}.
    - Rank 3 (non-planar): dim = 5 (1 component)
    - Rank 2 (planar, non-collinear): dim = 4 (2 components)
    - Rank 1 (collinear): dim = 2 (6 components)
    """

    # Dimension of the top-dimensional stratum (non-planar configurations)
    y1 = 5

    # Dimensions of the next strata (planar, non-collinear configurations)
    # This stratum has two connected components.
    y2 = 4
    y3 = 4

    # Dimensions of the lowest-dimensional strata (collinear configurations)
    # This stratum has 6 connected components.
    y_collinear = 2
    num_collinear_components = 6

    # Assemble the list of dimensions in descending order
    dimensions = [y1, y2, y3] + [y_collinear] * num_collinear_components

    # Print the result in the specified format
    print(",".join(map(str, dimensions)))

solve_configuration_space_dimensions()