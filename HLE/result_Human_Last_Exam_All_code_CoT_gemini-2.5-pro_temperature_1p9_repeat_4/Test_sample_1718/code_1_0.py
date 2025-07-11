def calculate_riemann_components_kahler():
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.
    """
    try:
        m_input = input("Enter the complex dimension of the Kähler manifold (m): ")
        m = int(m_input)

        if m <= 0:
            print("Error: The complex dimension must be a positive integer.")
            return

        # On a Kähler manifold of complex dimension m, the independent components
        # of the Riemann tensor correspond to the entries of a d x d Hermitian matrix,
        # where d = m(m+1)/2.
        # The number of independent real components in such a matrix is d^2.

        # Formula: N = ( (m * (m + 1)) / 2 )^2

        d = (m * (m + 1)) // 2
        num_components = d ** 2

        print(f"\nFor a Kähler manifold with complex dimension m = {m}:")
        print("The formula for the number of independent components is (m * (m + 1) / 2)^2.")
        print(f"First, we calculate d = m * (m + 1) / 2:")
        print(f"d = ({m} * ({m} + 1)) / 2 = {d}")
        print(f"Then, the number of components is d^2:")
        print(f"Result = {d}^2 = {num_components}")

    except ValueError:
        print("Error: Invalid input. Please enter an integer for the complex dimension.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    calculate_riemann_components_kahler()