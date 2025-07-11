def count_riemann_kahler_components(m: int):
    """
    Calculates and prints the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension 'm'.

    Args:
        m: The complex dimension of the Kähler manifold (must be a positive integer).
    """
    if not isinstance(m, int) or m <= 0:
        print(f"Error: Complex dimension 'm' must be a positive integer. Got: {m}")
        return None

    # The number of components is given by N^2, where N is the dimension of
    # the space of symmetric tensors of type (2,0), which is m*(m+1)/2.
    # We use integer division // to ensure the result is an integer.
    N = (m * (m + 1)) // 2
    result = N ** 2

    # The prompt requests that we output each number in the final equation.
    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {2*m}):")
    print(f"The formula for the number of independent components is (m * (m + 1) / 2)^2.")
    print(f"Step 1: Calculate the intermediate value N = (m * (m + 1)) / 2.")
    print(f"   N = ({m} * ({m} + 1)) / 2 = {m * (m+1)} / 2 = {N}")
    print(f"Step 2: Square N to get the final result.")
    print(f"   Total Components = N^2 = {N}^2 = {result}")
    print("-" * 40)
    return result

# --- Main Execution ---
# The number of components depends on the dimension.
# Let's calculate it for a few illustrative examples.
print("Calculating the number of independent Riemann tensor components for sample Kähler manifolds.")
print("======================================================================================\n")

# For m=1 (e.g., the Riemann sphere S^2), the result is 1, as for any 2D manifold.
count_riemann_kahler_components(1)

# For m=2 (e.g., a K3 surface), a common example in geometry and string theory.
count_riemann_kahler_components(2)

# For m=3 (e.g., a Calabi-Yau threefold), another important case.
count_riemann_kahler_components(3)