import math

def calculate_riemann_components_kahler():
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.
    """
    try:
        m_str = input("Enter the complex dimension of the Kähler manifold (m): ")
        m = int(m_str)
        if m < 0:
            print("Dimension cannot be negative.")
            return
        if m == 0:
            print("For a 0-dimensional manifold, there is 1 component (a single point).")
            return

        # The number of independent entries in a symmetric pair of indices from 1 to m.
        d = m * (m + 1) // 2

        # For a Kähler manifold, the Riemann tensor components form a d x d Hermitian matrix.
        # The number of independent real components in such a matrix is d^2.
        num_components = d * d

        print(f"\nFor a Kähler manifold of complex dimension m = {m} (real dimension n = {2*m}):")
        print("The number of independent components of the Riemann tensor is given by the formula (m*(m+1)/2)^2.")
        
        # Displaying the calculation step-by-step as requested
        print("\nCalculation:")
        print(f"Number of components = ({m} * ({m} + 1) / 2)^2")
        print(f"                     = ({m} * {m+1} / 2)^2")
        print(f"                     = ({m * (m+1)} / 2)^2")
        print(f"                     = ({d})^2")
        print(f"                     = {num_components}")

    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_riemann_components_kahler()