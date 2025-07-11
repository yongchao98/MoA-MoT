def solve():
    """
    Calculates and explains the minimum hidden-layer width for a shallow neural network
    to compute the squared norm of an N-dimensional vector.
    """
    # The target function is the squared norm: f(x) = ||x||^2 = sum(x_i^2 for i=1..N).
    # This can be seen as the sum of N independent functions g(x_i) = x_i^2.

    # The activation function is GeLU, which is monotonic. The target function g(z) = z^2
    # is not monotonic. To approximate a non-monotonic function by summing monotonic ones,
    # at least two are needed (one for the decreasing part, one for the increasing part).

    # For each dimension `i`, we can approximate the x_i^2 term by using a pair of neurons,
    # for example, one whose activation is a function of `x_i` and another whose
    # activation is a function of `-x_i`. This provides the necessary components to
    # build a parabolic shape.
    neurons_per_dimension = 2

    # Since there are N dimensions and the function is a sum of N independent parabolas,
    # we need to apply this principle N times.
    # Thus, the total minimum number of hidden neurons is 2 * N.

    final_coefficient = neurons_per_dimension
    final_variable = 'N'

    print("The minimum required hidden-layer width (H) to compute the squared norm of an N-dimensional vector is given by the formula:")
    # The final equation expresses the width H in terms of N.
    print(f"H = {final_coefficient} * {final_variable}")

    print("\nTo follow the instruction 'output each number in the final equation!':")
    # The only explicit number in the final equation "2 * N" is 2.
    print(f"The coefficient is: {final_coefficient}")

solve()