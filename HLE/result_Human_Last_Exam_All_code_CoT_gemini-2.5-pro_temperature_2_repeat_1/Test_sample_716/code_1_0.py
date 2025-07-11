import math

def solve_hausdorff_dimension():
    """
    Determines the Hausdorff dimension of the curve parametrized by:
    x(t) = sin(pi*t), y(t) = sin(t), z(t) = cos(2t).
    """

    print("Step 1: The curve is a 1-dimensional object as it is parametrized by a single variable, t.")
    print("x(t) = sin(pi*t)")
    print("y(t) = sin(t)")
    print("z(t) = cos(2t)\n")

    print("Step 2: We analyze the smoothness of the curve. The Hausdorff dimension of a set is equal")
    print("to its topological dimension if the set is sufficiently smooth (at least continuously differentiable, C^1).\n")

    print("Step 3: The component functions x(t), y(t), and z(t) are defined using sine and cosine,")
    print("which are infinitely differentiable functions. This means the curve is a smooth (C-infinity) 1-dimensional manifold.\n")

    print("Step 4: According to a standard theorem in fractal geometry, the Hausdorff dimension")
    print("of any continuously differentiable (C^1) m-dimensional manifold is m.\n")
    
    # The dimension of our manifold (a curve) is m=1.
    hausdorff_dimension = 1
    
    print("Step 5: For our curve, the dimension m is 1. Therefore, its Hausdorff dimension is also 1.")
    print(f"Final equation: Hausdorff Dimension = {hausdorff_dimension}")

if __name__ == '__main__':
    solve_hausdorff_dimension()
