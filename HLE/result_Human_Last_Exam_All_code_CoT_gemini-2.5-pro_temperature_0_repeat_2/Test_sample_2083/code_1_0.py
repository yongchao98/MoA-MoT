def solve_for_width():
    """
    Calculates the minimum hidden-layer width for a shallow neural network
    to compute the squared norm of an N-dimensional vector.
    """
    # N represents the dimension of the input vector.
    # The problem asks for the minimum hidden-layer width in terms of N.
    # Based on the analysis, the minimum width is 2*N.

    # Let's demonstrate the calculation for a sample value of N.
    # You can change this value to see the result for a different dimension.
    N = 10

    # The minimum width is calculated as 2 * N.
    result = 2 * N

    print("The theoretical minimum hidden-layer width to compute the squared norm of an N-dimensional vector is 2*N.")
    print(f"For an example where N = {N}, the minimum width is calculated as follows:")
    
    # The prompt asks to output each number in the final equation.
    print(f"2 * {N} = {result}")

solve_for_width()