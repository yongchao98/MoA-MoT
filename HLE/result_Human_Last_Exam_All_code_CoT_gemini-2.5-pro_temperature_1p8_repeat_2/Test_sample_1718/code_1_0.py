def calculate_kahler_riemann_components(m_values):
    """
    Calculates the number of independent components of the Riemann tensor
    for a Kähler manifold of given complex dimensions.

    Args:
        m_values (list): A list of integer complex dimensions to calculate for.
    """
    print("The number of independent real components of the Riemann tensor on a Kähler manifold")
    print("of complex dimension 'm' is given by the formula: N = (m * (m + 1) / 2)^2.\n")

    for m in m_values:
        if not isinstance(m, int) or m <= 0:
            print(f"Skipping invalid dimension: {m}. Dimension must be a positive integer.")
            continue
            
        print(f"--- For a complex dimension m = {m} ---")
        
        # Real dimension n
        n = 2 * m
        
        # The complex dimension of the space of symmetric (2,0)-tensors
        d_num = m * (m + 1)
        d = d_num // 2
        
        # The number of independent real components is d^2
        num_components = d * d
        
        # For comparison, the number of components in a general Riemannian manifold of same real dimension
        # n = 2m
        if n >= 2:
            general_components = (n**2 * (n**2 - 1)) // 12
            print(f"A general Riemannian manifold of the same real dimension n = {n} would have {general_components} components.")
        
        print("Step 1: Calculate the dimension 'd' of the space of symmetric (2,0)-tensors.")
        print(f"d = (m * (m + 1)) / 2 = ({m} * ({m} + 1)) / 2 = {d_num} / 2 = {d}")
        
        print("\nStep 2: The total number of components 'N' is the square of 'd'.")
        print(f"N = d^2 = {d}^2 = {num_components}\n")

# Example usage for complex dimensions 1, 2, and 3
# A complex dimension of 1 corresponds to a Riemann surface (real dim 2).
# A complex dimension of 2 corresponds to surfaces like K3 or CP^2 (real dim 4).
example_dimensions = [1, 2, 3, 4]
calculate_kahler_riemann_components(example_dimensions)