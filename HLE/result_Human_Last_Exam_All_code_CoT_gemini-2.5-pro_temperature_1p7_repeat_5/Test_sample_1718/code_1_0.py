def calculate_riemann_components_kahler():
    """
    Calculates the number of independent components of the Riemann tensor
    for a Kähler manifold of a given complex dimension 'm'.
    """
    try:
        m_str = input("Enter the complex dimension of the Kähler manifold (m): ")
        m = int(m_str)

        if m <= 0:
            print("\nError: The complex dimension must be a positive integer.")
            return

        # Real dimension n = 2m
        n = 2 * m

        # First, calculate N, the number of components of a symmetric tensor of rank 2 in m dimensions.
        # This is m*(m+1)/2.
        N = (m * (m + 1)) // 2

        # The number of independent real components is N^2.
        result = N * N

        print(f"\nFor a Kähler manifold with complex dimension m = {m}:")
        print(f"The real dimension is n = 2 * m = {n}.")

        print("\nThe calculation involves two steps:")
        print("1. Find the number of symmetric pairs of indices, N.")
        print("2. Square N to get the total number of independent components.")

        print("\nStep 1: Calculate N")
        print(f"N = (m * (m + 1)) / 2")
        print(f"N = ({m} * ({m} + 1)) / 2 = ({m} * {m+1}) / 2 = {N}")

        print("\nStep 2: Calculate the final result from N")
        print(f"Number of components = N^2")
        print(f"Number of components = {N}^2 = {result}")

    except ValueError:
        print("\nError: Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"\nAn error occurred: {e}")

if __name__ == '__main__':
    calculate_riemann_components_kahler()