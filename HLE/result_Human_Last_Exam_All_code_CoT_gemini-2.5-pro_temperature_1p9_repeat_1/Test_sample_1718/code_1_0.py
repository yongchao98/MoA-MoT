import math

def calculate_riemann_components(m):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    # The real dimension is n = 2m
    n = 2 * m
    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {n}):")

    # The number of components is determined by the dimension N of the space S^2(C^m).
    # N = m * (m + 1) / 2
    N = m * (m + 1) // 2

    # The total number of independent real components is N^2.
    result = N ** 2

    # Print the equation with all numbers
    print(f"Number of components = (m * (m + 1) / 2)^2 = ({m} * ({m} + 1) / 2)^2 = {N}^2 = {result}")

if __name__ == '__main__':
    # As the user did not specify a dimension, we will calculate for a common example, m=2.
    # A K3 surface is a well-known example of a Kähler manifold with m=2.
    complex_dimension = 2
    calculate_riemann_components(complex_dimension)
