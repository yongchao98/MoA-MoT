def solve_robot_arm_configuration():
    """
    Solves for the dimensions of the manifold components of the robot arm configuration space X_4.

    The configuration space X_4 of 4-segment unit-length robot arms that form a closed loop
    can be decomposed into a disjoint union of connected manifolds. This decomposition is
    based on the geometric degeneracies of the configurations.

    The dimensions of these component manifolds are determined as follows:
    1.  Generic configurations (non-planar) form a single connected manifold of dimension 5.
    2.  Planar, non-collinear configurations form three distinct connected manifolds. These
        are configurations where a pair of non-adjacent vectors sum to zero (forming
        parallelograms or antiparallelograms). Each of these manifolds has dimension 4.
    3.  Collinear configurations form three distinct connected manifolds. These occur when
        two vectors are u and two are -u for some unit vector u. Each of these manifolds
        is a 2-sphere and thus has dimension 2.

    This gives a total of 1 + 3 + 3 = 7 disjoint connected manifolds.
    """

    # The dimensions of the disjoint connected manifolds, sorted in descending order.
    # y_1: The dimension of the main component of generic (non-planar) configurations.
    # y_2, y_3, y_4: The dimensions of the three components of planar, non-collinear configurations.
    # y_5, y_6, y_7: The dimensions of the three components of collinear configurations.
    y = [5, 4, 4, 4, 2, 2, 2]

    # The problem asks for the tuple (y_1, ..., y_l)
    # We print the numbers separated by commas.
    print(','.join(map(str, y)))

solve_robot_arm_configuration()