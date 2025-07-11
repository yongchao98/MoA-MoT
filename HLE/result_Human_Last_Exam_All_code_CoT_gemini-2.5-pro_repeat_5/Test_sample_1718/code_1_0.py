def count_kahler_riemann_components(m):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    print(f"Calculating for a Kähler manifold of complex dimension m = {m}")
    print(f"The real dimension is n = 2 * m = {2 * m}\n")

    # Step 1: Calculate the number of symmetric pairs of indices.
    # The formula is N_sym = m * (m + 1) / 2
    print("Step 1: Calculate the number of symmetric pairs of indices (N_sym).")
    print(f"N_sym = m * (m + 1) / 2")
    m_plus_1 = m + 1
    numerator = m * m_plus_1
    n_sym = numerator // 2
    print(f"N_sym = {m} * ({m} + 1) / 2 = {m} * {m_plus_1} / 2 = {numerator} / 2 = {n_sym}\n")

    # Step 2: The total number of independent components is N_sym^2.
    print("Step 2: The total number of independent components is N_sym squared.")
    print("Result = N_sym^2")
    result = n_sym ** 2
    print(f"Result = {n_sym}^2 = {result}\n")

    print(f"A Kähler manifold of complex dimension {m} has {result} independent Riemann tensor components.")

# --- Main execution ---
# We will calculate the number of components for a complex dimension m = 2.
# This corresponds to a real 4-dimensional manifold, like a K3 surface.
if __name__ == "__main__":
    complex_dimension = 2
    count_kahler_riemann_components(complex_dimension)
