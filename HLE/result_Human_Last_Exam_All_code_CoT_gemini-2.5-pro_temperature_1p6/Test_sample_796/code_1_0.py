def solve_robot_arm_dimensions():
    """
    Calculates the dimensions of the manifolds in the decomposition of the
    configuration space X_4 of a 4-segment robot arm.

    The configuration space X_4 consists of 4-tuples of unit vectors
    (v1, v2, v3, v4) that sum to the zero vector. The space is decomposed
    into disjoint connected manifolds by stratifying based on the length
    L = ||v1 + v2||.
    """

    # Dimension of the building block of our space, the 2-sphere S^2.
    dim_S2 = 2

    # The total number of vectors in the arm.
    num_vectors = 4

    # --- Dimension of the generic manifold (y1) ---
    # This corresponds to the case where 0 < L < 2.
    # The dimension is the total degrees of freedom of the four vectors minus
    # the number of constraints imposed by the sum-to-zero equation.
    # Total DoF = 4 * dim(S^2) = 4 * 2 = 8.
    # The vector equation v1+v2+v3+v4 = 0 provides 3 independent constraints.
    # Therefore, the dimension y1 = 8 - 3 = 5.
    y1 = num_vectors * dim_S2 - 3
    print(f"Dimension of the generic manifold (0 < L < 2): {y1}")


    # --- Dimension of the first degenerate manifold (y2) ---
    # This corresponds to the case where L = ||v1 + v2|| = 0.
    # This implies v2 = -v1. The governing equation becomes v3 + v4 = 0,
    # so v4 = -v3. The configuration is determined by an independent choice
    # of v1 in S^2 and v3 in S^2.
    # The space is homeomorphic to S^2 x S^2.
    # Its dimension y2 = dim(S^2) + dim(S^2) = 2 + 2 = 4.
    y2 = dim_S2 + dim_S2
    print(f"Dimension of the L=0 manifold (S^2 x S^2): {y2}")


    # --- Dimension of the second degenerate manifold (y3) ---
    # This corresponds to the case where L = ||v1 + v2|| = 2.
    # This implies v2 = v1. The governing equation becomes v3 + v4 = -2*v1.
    # The length ||v3+v4|| must be 2, which implies v3 = v4.
    # The equation simplifies to 2*v3 = -2*v1, so v3 = -v1.
    # The configuration is fully determined by the choice of v1 in S^2.
    # The space is homeomorphic to S^2.
    # Its dimension y3 = dim(S^2) = 2.
    y3 = dim_S2
    print(f"Dimension of the L=2 manifold (S^2): {y3}")


    # --- Final Answer ---
    # The dimensions of the three disjoint manifolds are 5, 4, and 2.
    # We list them in descending order as requested.
    dimensions = sorted([y1, y2, y3], reverse=True)
    
    print("\nThe dimensions (y1, ..., yl) for the least l are:")
    # "Remember in the final code you still need to output each number in the final equation!"
    # The following line prints each number in the final result tuple.
    print(f"{dimensions[0]},{dimensions[1]},{dimensions[2]}")

solve_robot_arm_dimensions()