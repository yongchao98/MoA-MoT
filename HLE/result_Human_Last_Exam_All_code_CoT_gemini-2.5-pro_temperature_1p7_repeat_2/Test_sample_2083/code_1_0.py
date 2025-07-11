def solve_network_width():
    """
    Calculates the minimum hidden-layer width for a shallow neural network
    to compute the squared norm of an N-dimensional vector.
    """
    # Let N be the dimension of the input vector.
    # We will use an example value for N.
    N = 10
    
    # The minimum hidden-layer width H required is 2*N.
    # Lower bound: Approximating an even function with GeLU requires pairing neurons
    # to cancel odd terms. Analysis of the Hessian matrix shows that N pairs
    # of neurons are necessary. This gives H >= 2*N.
    # Upper bound: A constructive proof shows that a network with 2*N neurons
    # can approximate the squared norm to arbitrary precision.
    # The construction sums N approximations of x_i^2, where each approximation
    # uses 2 neurons: GeLU(s*x_i) and GeLU(-s*x_i).
    # Since the lower bound and upper bound match, the minimum width is 2*N.
    
    minimum_width = 2 * N
    
    print("For a shallow neural network with one hidden layer using the GeLU activation function,")
    print("the minimum hidden-layer width H required to compute the squared norm of an N-dimensional input vector is H = 2*N.")
    print("\nFor an example input dimension:")
    print(f"N = {N}")
    print("\nThe final equation is:")
    print(f"H = 2 * {N} = {minimum_width}")

solve_network_width()
