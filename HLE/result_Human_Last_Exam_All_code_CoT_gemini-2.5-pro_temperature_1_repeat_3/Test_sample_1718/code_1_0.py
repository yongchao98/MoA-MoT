def calculate_kahler_riemann_components():
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold.
    """
    try:
        n_str = input("Enter the real dimension of the Kähler manifold (must be an even positive integer): ")
        n = int(n_str)
        if n <= 0 or n % 2 != 0:
            print("Error: The real dimension of a Kähler manifold must be a positive even integer.")
            return
    except ValueError:
        print("Error: Please enter a valid integer.")
        return

    # Complex dimension
    m = n // 2

    # Number of ways to choose a symmetric pair of indices from m dimensions
    intermediate_numerator = m * (m + 1)
    intermediate_base = intermediate_numerator // 2
    
    # Final result is the square of the above
    result = intermediate_base ** 2

    # Print the explanation and step-by-step calculation
    print(f"\nFor a Kähler manifold of real dimension n = {n}:")
    print(f"The complex dimension is m = n / 2 = {m}.")
    print("The number of independent components is given by the formula: (m * (m + 1) / 2)^2")
    print("Calculation:")
    print(f"(({m} * ({m} + 1)) / 2)^2 = (({m} * {m + 1}) / 2)^2 = ({intermediate_numerator} / 2)^2 = {intermediate_base}^2 = {result}")

if __name__ == '__main__':
    calculate_kahler_riemann_components()