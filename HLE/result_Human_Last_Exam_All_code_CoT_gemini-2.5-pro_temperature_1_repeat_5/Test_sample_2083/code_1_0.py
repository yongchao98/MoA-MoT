def solve_network_width():
    """
    Calculates and explains the minimum hidden-layer width for a shallow
    neural network with GeLU activation to compute the squared norm of an
    N-dimensional vector.
    """

    # The problem is to find the minimum hidden-layer width (H) to compute
    # f(x) = ||x||^2 = x_1^2 + x_2^2 + ... + x_N^2.

    # 1. The function f(x) is a sum of N one-dimensional quadratic functions.
    #    We can dedicate a set of neurons to compute each x_i^2 term.

    # 2. To compute z^2, we need to analyze the GeLU activation function.
    #    The function z^2 is an even function. GeLU is not.
    #    Approximating an even function requires at least two GeLU neurons
    #    to create a symmetric response. For example, c*(GeLU(z) + GeLU(-z)).
    neurons_per_dimension = 2

    # 3. N is the number of dimensions of the input vector.
    #    The final answer is expressed in terms of N.
    num_dimensions_symbol = "N"

    # 4. The total minimum width H is the number of neurons per dimension
    #    multiplied by the number of dimensions.
    #    H = 2 * N

    print("To compute the squared norm of an N-dimensional vector, the network must approximate the function:")
    print("f(x) = x_1^2 + x_2^2 + ... + x_N^2\n")

    print(f"A single term, z^2, can be approximated using {neurons_per_dimension} GeLU neurons.")
    print("This is because z^2 is an even function, while a single GeLU neuron is not, so a combination of at least two is required.\n")

    print("To compute the sum for an N-dimensional vector, we need to apply this principle for each dimension.")
    print("This results in a total minimum number of hidden neurons (H) given by the equation:\n")

    # Output the final equation with each component, as requested.
    coefficient = neurons_per_dimension
    variable = num_dimensions_symbol
    print(f"H_min = {coefficient} * {variable}")

solve_network_width()
<<<2*N>>>