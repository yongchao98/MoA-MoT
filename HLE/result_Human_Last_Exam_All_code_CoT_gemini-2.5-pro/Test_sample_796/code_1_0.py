def solve_robot_arm_dimensions():
    """
    This function determines the dimensions of the manifold components of the robot arm configuration space.
    The space X_4 of 4-segment unit-length arms forming a closed loop is decomposed into a disjoint union
    of connected manifolds without boundary. This decomposition is based on the geometric properties of the configurations.

    The decomposition yields:
    - One component of dimension 5, consisting of generic non-planar configurations.
    - Three distinct components of dimension 4, consisting of configurations where one pair of vectors is antipodal, but the configuration is not fully collinear.
    - Three distinct components of dimension 2, consisting of the fully collinear configurations.

    This results in a total of 1 + 3 + 3 = 7 manifold components.
    The dimensions are y_1=5, y_2=4, y_3=4, y_4=4, y_5=2, y_6=2, y_7=2.
    """
    
    # The dimensions of the manifold components, sorted in descending order.
    dimensions = [5, 4, 4, 4, 2, 2, 2]
    
    # The problem asks to output each number in the final equation.
    # We will print them in the specified format.
    print(','.join(map(str, dimensions)))

solve_robot_arm_dimensions()