def calculate_kahler_riemann_components():
    """
    Calculates the number of independent real components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.

    The formula is N^2, where N = m * (m + 1) / 2 and 'm' is the complex dimension.
    """
    try:
        # Prompt the user to enter the complex dimension 'm'.
        m_str = input("Enter the complex dimension 'm' of the Kähler manifold: ")
        m = int(m_str)

        if m <= 0:
            print("Error: The complex dimension must be a positive integer.")
            return

        # Step 1: Calculate N, the number of symmetric index pairs.
        # This is the dimension of the space of symmetric m x m matrices.
        N = m * (m + 1) // 2

        # Step 2: The number of independent components is N^2, as the
        # relevant components form an N x N Hermitian matrix.
        num_components = N * N

        # Print the results, showing the calculation as requested.
        print(f"\nFor a Kähler manifold of complex dimension m = {m}:")
        print(f"The number of symmetric index pairs is N = (m * (m + 1)) / 2 = ({m} * ({m} + 1)) / 2 = {N}")
        print(f"The number of independent real components of the Riemann tensor is N^2 = {N}^2 = {num_components}")

    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    calculate_kahler_riemann_components()