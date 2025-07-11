def calculate_kahler_riemann_components(m: int):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    Args:
        m: The complex dimension of the Kähler manifold.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    # The formula is (m * (m + 1) / 2)^2
    # We will show the calculation step-by-step as requested.
    
    print(f"--- Calculating for a Kähler manifold of complex dimension m = {m} ---")
    
    # Step 1: Calculate m + 1
    m_plus_1 = m + 1
    print(f"The first step in the numerator is m + 1: {m} + 1 = {m_plus_1}")
    
    # Step 2: Calculate the full numerator m * (m + 1)
    numerator = m * m_plus_1
    print(f"The full numerator is m * (m + 1): {m} * {m_plus_1} = {numerator}")

    # Step 3: Calculate the base of the exponent
    base = numerator / 2
    # Ensure the base is an integer for clarity in the printout
    base_int = int(base)
    print(f"The base of the exponent is (m * (m + 1)) / 2: {numerator} / 2 = {base_int}")

    # Step 4: Calculate the final result by squaring the base
    final_count = base_int ** 2
    print(f"The final number of independent components is the base squared: {base_int}^2 = {final_count}")
    print("-" * (60))
    

if __name__ == '__main__':
    # The real dimension n of the manifold is 2*m.
    
    # Example 1: A Riemann surface (like a torus or a sphere). Complex dimension m=1, real dimension n=2.
    calculate_kahler_riemann_components(1)

    # Example 2: A K3 surface or complex projective plane CP^2. Complex dimension m=2, real dimension n=4.
    calculate_kahler_riemann_components(2)

    # Example 3: A Calabi-Yau threefold (important in string theory). Complex dimension m=3, real dimension n=6.
    calculate_kahler_riemann_components(3)
    
    # Example 4: Complex dimension m=4, real dimension n=8
    calculate_kahler_riemann_components(4)