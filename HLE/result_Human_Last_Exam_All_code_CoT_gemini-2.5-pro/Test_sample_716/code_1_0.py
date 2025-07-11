def find_hausdorff_dimension():
    """
    This script determines the Hausdorff dimension of the curve parametrized by:
    x(t) = sin(pi*t)
    y(t) = sin(t)
    z(t) = cos(2*t)

    The solution is based on a mathematical theorem about smooth curves.
    """

    # Step 1: Understand the principle.
    # The Hausdorff dimension of a C^1-smooth curve (a curve with a continuous
    # non-zero derivative) is equal to its topological dimension, which is 1.
    # We need to verify if the given curve is C^1-smooth.

    # Step 2: Analyze the curve's parametrization.
    # The parametric functions x(t), y(t), and z(t) are combinations of sine and
    # cosine functions, which are infinitely differentiable (and thus C^1).

    # Step 3: Analyze the tangent vector r'(t) = (dx/dt, dy/dt, dz/dt).
    # The derivatives are:
    # dx/dt = pi * cos(pi*t)
    # dy/dt = cos(t)
    # dz/dt = -2 * sin(2*t)
    # These derivatives are all continuous functions.

    # Step 4: Check if the tangent vector can be the zero vector.
    # This requires all three components to be zero for the same value of t.
    # 1. pi * cos(pi*t) = 0  =>  pi*t = pi/2 + k*pi  =>  t = k + 1/2 for integer k
    # 2. cos(t) = 0         =>  t = pi/2 + m*pi for integer m (t approx 1.57, 4.71, ...)
    # 3. -2 * sin(2*t) = 0      =>  2*t = n*pi       =>  t = n*pi/2 for integer n

    # Let's check if these conditions can hold simultaneously.
    # From condition 1, t must be a half-integer (e.g., 0.5, 1.5, -0.5).
    # If we take t = 0.5 from condition 1, condition 2 gives cos(0.5), which is not 0.
    # If we take a value for t from condition 2, like t = pi/2, then condition 1 gives
    # pi * cos(pi^2 / 2), which is not 0.
    # Therefore, there is no value of 't' for which the tangent vector is (0, 0, 0).
    # The curve is C^1-smooth.

    # Step 5: State the conclusion.
    # Since the curve is C^1-smooth, its Hausdorff dimension is 1.
    hausdorff_dimension = 1

    print("The analysis shows that the curve is C^1-smooth.")
    print("For any C^1-smooth curve, the Hausdorff dimension is equal to its topological dimension.")
    print("The final equation for the Hausdorff Dimension (D_H) is:")
    # The final equation has one number: 1
    print(f"D_H = {hausdorff_dimension}")


if __name__ == "__main__":
    find_hausdorff_dimension()