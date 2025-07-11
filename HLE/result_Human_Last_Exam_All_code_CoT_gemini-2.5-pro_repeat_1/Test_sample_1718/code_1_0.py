def calculate_kahler_riemann_components(m):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    # The real dimension n is 2*m
    n = 2 * m

    # On a Kähler manifold, the Riemann tensor is determined by a Hermitian form on the
    # space of symmetric (2,0)-tensors. The dimension of this space is N = C(m+1, 2).
    # N = m * (m + 1) / 2
    N = m * (m + 1) // 2

    # The number of independent real components in an N x N Hermitian matrix is N^2.
    num_components = N ** 2

    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {n}):")
    print("\nStep 1: The number of components is given by the formula ( (m * (m + 1)) / 2 )^2.")

    print("\nStep 2: Substitute the value of m = {m} into the formula.")
    # Show the numbers in the final equation
    m_plus_1 = m + 1
    numerator = m * m_plus_1
    
    print(f"Equation: ( ({m} * ({m} + 1)) / 2 )^2")
    print(f"        = ( ({m} * {m_plus_1}) / 2 )^2")
    print(f"        = ( {numerator} / 2 )^2")
    print(f"        = ( {N} )^2")
    print(f"        = {num_components}")
    
    print(f"\nResult: There are {num_components} independent entries in the Riemann tensor.")


# Example calculation for a Kähler manifold of complex dimension m = 2 (e.g., a K3 surface)
complex_dimension_m = 2
calculate_kahler_riemann_components(complex_dimension_m)