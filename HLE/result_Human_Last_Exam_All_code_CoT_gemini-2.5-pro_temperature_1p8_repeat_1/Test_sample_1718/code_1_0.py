def count_riemann_components_kahler(m):
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of complex dimension m.

    Args:
        m (int): The complex dimension of the Kähler manifold.
    
    Returns:
        int: The number of independent components.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: The complex dimension 'm' must be a positive integer.")
        return None

    # The formula is (m * (m + 1) / 2)^2.
    # The prompt asks to output each number in the final equation.
    # We will show the calculation step-by-step.
    
    term1 = m
    term2 = m + 1
    divisor = 2
    
    # Calculate the base of the exponent, N
    base = (term1 * term2) // divisor
    
    # Square the result to get the final count
    num_components = base ** 2

    print(f"For a Kähler manifold of complex dimension m = {m} (real dimension n = {2*m}):")
    print("The number of independent components is given by the formula: (m * (m + 1) / 2)^2")
    
    # Printing each number in the equation for clarity
    print("\nCalculation steps:")
    print(f"1. Calculate N = (m * (m + 1)) / 2")
    print(f"   N = ({term1} * ({term1} + 1)) / {divisor}")
    print(f"   N = ({term1} * {term2}) / {divisor}")
    print(f"   N = {term1 * term2} / {divisor}")
    print(f"   N = {base}")
    
    print(f"\n2. Calculate the final number = N^2")
    print(f"   Number of components = {base}^2")
    print(f"   Number of components = {num_components}")
    print("-" * 40)
    
    return num_components

if __name__ == "__main__":
    print("Calculating the number of independent Riemann tensor components for example Kähler manifolds.\n")
    
    # Case m=1 (e.g., a Riemann surface, real dimension 2)
    count_riemann_components_kahler(1)
    
    # Case m=2 (e.g., a K3 surface, real dimension 4)
    count_riemann_components_kahler(2)
    
    # Case m=3 (e.g., a Calabi-Yau threefold, real dimension 6)
    count_riemann_components_kahler(3)