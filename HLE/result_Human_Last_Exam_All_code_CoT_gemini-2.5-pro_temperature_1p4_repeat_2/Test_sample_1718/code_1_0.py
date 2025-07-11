def calculate_riemann_components_kahler(m):
    """
    Calculates and prints the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension 'm'.

    The function prints the calculation process step-by-step, including the final equation
    with all numbers substituted.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return

    print(f"For a Kähler manifold of complex dimension m = {m}:")
    
    # The formula is N = (k)^2, where k = (m * (m + 1)) / 2.
    
    # Step 1: Calculate the term m * (m + 1)
    m_plus_1 = m + 1
    k_numerator = m * m_plus_1
    
    # Step 2: Calculate k. Since m or m+1 is even, this is always an integer.
    k = k_numerator // 2
    
    # Step 3: Calculate the final result, N = k^2.
    result = k ** 2
    
    # Print the equation with all numbers substituted, following the calculation steps.
    print(f"The number of independent components N is given by the formula: ((m * (m + 1)) / 2)^2")
    print(f"Substituting m = {m}:")
    print(f"N = (({m} * ({m} + 1)) / 2)^2")
    print(f"N = (({m} * {m_plus_1}) / 2)^2")
    print(f"N = ({k_numerator} / 2)^2")
    print(f"N = {k}^2")
    print(f"N = {result}")
    print("-" * 30)

# We can now use this function to find the number of components for a few example dimensions.
# For m=1 (a Riemann surface, real dimension n=2)
calculate_riemann_components_kahler(1)

# For m=2 (a K3 surface, for example, real dimension n=4)
calculate_riemann_components_kahler(2)

# For m=3 (a Calabi-Yau manifold, for example, real dimension n=6)
calculate_riemann_components_kahler(3)
