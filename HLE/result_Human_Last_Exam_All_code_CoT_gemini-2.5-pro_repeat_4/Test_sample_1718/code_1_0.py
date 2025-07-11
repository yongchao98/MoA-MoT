def count_riemann_tensor_kahler(m: int):
    """
    Calculates the number of independent entries of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m: The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: Please provide a positive integer for the complex dimension 'm'.")
        return

    print(f"For a Kähler manifold of complex dimension m = {m}:")
    print("-" * 30)
    print("Step 1: Find the number of symmetric pairs of indices, N.")
    print(f"The formula is N = (m * (m + 1)) / 2")
    n_pairs_val = (m * (m + 1)) // 2
    print(f"N = ({m} * ({m} + 1)) / 2 = ({m} * {m+1}) / 2 = {m * (m+1)} / 2 = {n_pairs_val}")
    print("\nStep 2: The number of independent components is N^2.")
    print(f"Number of components = N^2")
    result = n_pairs_val ** 2
    print(f"Number of components = {n_pairs_val}^2 = {result}")
    print("-" * 30)
    print(f"Final Equation: ( ({m} * ({m} + 1)) / 2 )^2 = {result}")

# Example: Calculate for a manifold of complex dimension m=2 (like a K3 surface)
# You can change this value to any positive integer.
complex_dimension = 2
count_riemann_tensor_kahler(complex_dimension)