def solve_task():
    """
    Calculates the minimum hidden-layer width required for a shallow neural network
    with GeLU activations to compute the squared norm of an N-dimensional vector.
    """
    # Let N be the dimension of the input vector.
    # For demonstration, we'll use N = 10.
    N = 10

    # The analysis shows that the minimum number of hidden neurons (H) required
    # is given by the formula: H = N + 1.
    # This is derived by matching the network's properties (value, gradient, Hessian)
    # to the target function ||x||^2 at the origin, which leads to a set of
    # constraints on the network's weights and biases. The minimum number of
    # neurons needed to satisfy these constraints is N + 1.

    # The constant term in the equation is 1.
    constant_term = 1

    # Calculate the minimum width.
    minimum_width = N + constant_term

    # As requested, we output the numbers from the final equation: H = N + 1
    print(f"For an N-dimensional input vector, where N = {N},")
    print(f"the minimum hidden-layer width (H) is calculated by the equation:")
    print(f"H = {N} + {constant_term}")
    print(f"Result: The minimum width required is {minimum_width}.")


solve_task()