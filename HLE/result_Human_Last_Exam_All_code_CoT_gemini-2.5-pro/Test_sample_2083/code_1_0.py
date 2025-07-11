import sys

def calculate_minimum_width(N):
    """
    Calculates the minimum hidden-layer width (H) required for a shallow
    neural network to compute the squared norm of an N-dimensional input vector.

    The problem is to find the minimum H such that a network with the architecture
    N -> H -> 1 (with GeLU activations) can approximate f(x) = ||x||^2.

    The derivation shows that to synthesize the required quadratic form while satisfying
    the necessary mathematical constraints, the number of hidden neurons H must be
    at least N + 1.

    Args:
        N (int): The dimension of the input vector.

    Returns:
        None. Prints the result of the calculation.
    """
    if not isinstance(N, int) or N <= 0:
        print("Error: The input dimension N must be a positive integer.", file=sys.stderr)
        return

    # The minimum hidden-layer width H is N + 1.
    H = N + 1

    print(f"For an input vector of dimension N = {N}, the minimum required hidden-layer width is H.")
    print("The formula derived from the network's approximation requirements is:")
    # Output each number in the final equation as requested.
    print(f"H = N + 1 = {N} + 1 = {H}")

if __name__ == '__main__':
    # You can change this value to test with different input dimensions.
    input_dimension_N = 10
    calculate_minimum_width(input_dimension_N)