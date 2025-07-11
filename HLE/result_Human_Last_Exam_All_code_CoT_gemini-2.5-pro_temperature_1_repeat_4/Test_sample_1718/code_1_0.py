def count_kahler_riemann_components(m):
    """
    Calculates the number of independent real components of the Riemann tensor
    on a Kähler manifold of complex dimension m.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    print(f"Calculating for a Kähler manifold of complex dimension m = {m}")
    print("-" * 50)
    
    # Step 1: Calculate the number of symmetric pairs of indices.
    # This corresponds to the dimension of the space of symmetric m x m matrices.
    # N = m * (m + 1) / 2
    N_numerator = m * (m + 1)
    N = N_numerator // 2
    
    print("The number of independent components is determined by the number of symmetric pairs of holomorphic indices, N.")
    print(f"The formula for N is: m * (m + 1) / 2")
    print(f"N = {m} * ({m} + 1) / 2 = {m} * {m+1} / 2 = {N}")

    # Step 2: The total number of independent real components is N^2.
    # This is because the components form an N x N Hermitian matrix.
    num_components = N ** 2
    
    print("\nThe total number of independent real components is N^2.")
    print(f"Number of components = {N}^2 = {num_components}")
    print("-" * 50)
    print()

# --- Main execution ---
# Demonstrate the calculation for a few example dimensions.

# For m=1 (a Riemann surface like the 2-sphere)
# The real dimension is n=2.
count_kahler_riemann_components(1)

# For m=2 (e.g., a K3 surface)
# The real dimension is n=4.
count_kahler_riemann_components(2)

# For m=3 (e.g., a Calabi-Yau threefold)
# The real dimension is n=6.
count_kahler_riemann_components(3)
