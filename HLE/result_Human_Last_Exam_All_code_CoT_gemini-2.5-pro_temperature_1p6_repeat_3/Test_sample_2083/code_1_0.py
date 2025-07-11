def solve():
    """
    Calculates and prints the formula for the minimum hidden-layer width (H)
    of a shallow neural network required to compute the squared norm of an
    N-dimensional input vector.
    """

    # N represents the dimensionality of the input vector. It's a variable.
    variable_name = "N"

    # The coefficient is derived from the number of GELU neurons needed to
    # approximate a single quadratic term (z^2).
    # A single GELU neuron is not symmetric and cannot approximate z^2.
    # A pair of GELU neurons can be combined to form a symmetric function that
    # provides a local quadratic approximation. Thus, 2 neurons are needed for each
    # of the N dimensions.
    coefficient = 2

    print("To compute the squared norm of an N-dimensional vector (sum of x_i^2), the network must approximate N independent quadratic functions.")
    print("Approximating a single quadratic function x_i^2 requires a minimum of 2 GELU neurons.")
    print("Therefore, for N dimensions, the minimum required hidden-layer width (H) is:")
    print(f"H = {coefficient} * {variable_name}")
    print("\n--- Equation Breakdown ---")
    print(f"The formula for the minimum width H is a product of two numbers:")
    print(f"1. The coefficient: {coefficient}")
    print(f"2. The input dimension (a variable): {variable_name}")

solve()