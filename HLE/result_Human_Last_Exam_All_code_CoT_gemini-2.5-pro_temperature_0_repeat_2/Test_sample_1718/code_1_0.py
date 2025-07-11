def calculate_kahler_riemann_components(m):
    """
    Calculates the number of independent real components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    print(f"The formula for the number of independent components of the Riemann tensor on a Kähler manifold of complex dimension 'm' is (m * (m + 1) / 2)^2.")
    print(f"Let's calculate this for m = {m}.")
    print("-" * 30)

    # Step 1: Calculate the term inside the parenthesis
    # This term, d = m*(m+1)/2, is the dimension of the space of symmetric m x m matrices.
    d_numerator = m * (m + 1)
    d = d_numerator // 2

    # Step 2: Square the result
    # The total number of components is d^2.
    result = d * d

    # Print the final equation with all numbers, as requested
    print("The calculation is as follows:")
    print(f"({m} * ({m} + 1) / 2)^2 = ({m} * {m + 1} / 2)^2 = ({d_numerator} / 2)^2 = ({d})^2 = {result}")

    print("-" * 30)
    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {2*m}), the Riemann tensor has {result} independent entries.")

# --- Main execution ---
# We will use a common example: a manifold of complex dimension m = 2.
# This corresponds to a real 4-dimensional manifold.
# For comparison, a general 4D Riemannian manifold has 20 independent components.
complex_dimension_m = 2
calculate_kahler_riemann_components(complex_dimension_m)