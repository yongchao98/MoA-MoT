def solve_robot_arm_dimensions():
    """
    This function provides the dimensions of the disjoint connected manifolds
    that constitute the configuration space X_4.

    The configuration space X_4 is for a 4-segment unit-length robot arm
    that starts and ends at the origin. This space is known to be a singular
    algebraic variety. The problem asks for its decomposition into a disjoint
    union of connected manifolds, which corresponds to its Whitney stratification.

    The dimensions are determined through mathematical analysis:
    - y_1 = 5: This is the dimension of the main stratum, consisting of
      non-planar quadrilaterals. It is calculated as dim((S^2)^4) - 3 = 8 - 3 = 5.
    - y_2 = 4: This is the dimension of the first singular stratum, consisting
      of planar (but not collinear) quadrilaterals. This adds one constraint,
      reducing the dimension by 1.
    - y_3, y_4, y_5 = 2: This is the dimension of the most singular stratum,
      where all vectors are collinear. This stratum is not connected and splits
      into 3 disjoint components (homeomorphic to S^2), corresponding to the three
      ways to pair up opposing vectors: (v1,v2)(v3,v4), (v1,v3)(v2,v4),
      and (v1,v4)(v2,v3).
    """

    # The dimensions of the l=5 pairwise disjoint connected manifolds.
    y1 = 5
    y2 = 4
    y3 = 2
    y4 = 2
    y5 = 2

    # The final answer requires printing each number in the sequence.
    print(f"{y1},{y2},{y3},{y4},{y5}")

solve_robot_arm_dimensions()