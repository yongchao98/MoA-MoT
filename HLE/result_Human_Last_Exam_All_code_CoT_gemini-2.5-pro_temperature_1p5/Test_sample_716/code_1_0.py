def get_hausdorff_dimension():
    """
    Determines and explains the Hausdorff dimension of the given curve.

    The curve is defined in Cartesian coordinates by:
    x(t) = sin(pi * t)
    y(t) = sin(t)
    z(t) = cos(2t)
    """

    print("Step 1: Analyze the parametric functions.")
    print("The functions x(t) = sin(pi*t), y(t) = sin(t), and z(t) = cos(2t) are all infinitely differentiable (smooth).")
    print("-" * 20)

    print("Step 2: Relate smoothness to dimension.")
    print("A curve defined by smooth functions is a smooth 1-dimensional manifold, as long as it isn't just a point.")
    print("The tangent vector r'(t) = (pi*cos(pi*t), cos(t), -2*sin(2t)) is never the zero vector, which confirms the object is a proper curve.")
    print("-" * 20)

    print("Step 3: State the relevant mathematical principle.")
    print("For any smooth k-dimensional manifold, the Hausdorff dimension is equal to its topological dimension, which is k.")
    print("-" * 20)

    print("Step 4: Conclude the dimension.")
    print("Since our object is a smooth 1-dimensional curve (a 1-manifold), its Hausdorff dimension is 1.")
    print("-" * 20)

    # The final "equation" is the statement of the dimension's value.
    dimension = 1
    print(f"Final Result: The Hausdorff Dimension = {dimension}")

if __name__ == "__main__":
    get_hausdorff_dimension()