import math

def calculate_riemann_components_kahler(n):
    """
    Calculates the number of independent real components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        n (int): The complex dimension of the Kähler manifold.
                 The real dimension is 2n.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The complex dimension 'n' must be a positive integer.")
        return

    print(f"Calculating the number of independent Riemann tensor components for a Kähler manifold of complex dimension n = {n}.\n")

    # The formula is ((n * (n + 1)) / 2)^2
    # Step 1: Calculate the dimension of the space of symmetric tensors, N.
    print("Step 1: The components can be described by a Hermitian form on the space of symmetric n x n complex tensors.")
    print("         The dimension of this space, let's call it N, is calculated as:")
    print("N = (n * (n + 1)) / 2")
    n_plus_1 = n + 1
    numerator = n * n_plus_1
    N_float = numerator / 2
    N = int(N_float) # This division will always result in an integer

    print(f"N = ({n} * ({n} + 1)) / 2")
    print(f"N = ({n} * {n_plus_1}) / 2")
    print(f"N = {numerator} / 2")
    print(f"N = {N}\n")


    # Step 2: The number of independent real components is N^2.
    print("Step 2: The number of independent real components of a Hermitian form on an N-dimensional space is N^2.")
    print("Total Components = N^2")
    result = N**2
    print(f"Total Components = {N}^2")
    print(f"Total Components = {result}\n")

    print(f"Result: A Kähler manifold of complex dimension n={n} (real dimension {2*n}) has {result} independent Riemann tensor components.")

# --- User-configurable variable ---
# Set the complex dimension 'n' for the calculation.
# For example:
# n = 1 corresponds to a 2D Riemann surface (e.g., a sphere or torus).
# n = 2 corresponds to a 4D K3 surface or a complex projective plane.
# n = 3 corresponds to a 6D Calabi-Yau manifold.
complex_dimension_n = 3

calculate_riemann_components_kahler(complex_dimension_n)
