def solve_network_width():
    """
    This function determines the minimum hidden-layer width for a shallow neural network
    with GELU activation to compute the squared norm of an N-dimensional vector.
    """

    # As derived in the plan, the function ||x||^2 = Σ(x_i^2) is a sum of N independent even functions (z^2).
    # To approximate an even function z^2 using the non-even GELU activation, we need
    # a symmetric pairing of neurons. The simplest such pairing for each dimension x_i
    # is GELU(w*x_i) and GELU(-w*x_i), which requires 2 neurons.
    # Since there are N dimensions, we need 2 neurons for each.

    # The minimum hidden-layer width H is therefore 2 * N.
    # We will print the components of this equation.

    coefficient = 2
    variable = 'N'

    print("The minimum hidden-layer width, denoted by H, is determined by the dimensionality of the input vector, N.")
    print("The relationship is given by the following equation:")
    # Using f-string to clearly display the formula's components
    print(f"H = {coefficient} * {variable}")

solve_network_width()