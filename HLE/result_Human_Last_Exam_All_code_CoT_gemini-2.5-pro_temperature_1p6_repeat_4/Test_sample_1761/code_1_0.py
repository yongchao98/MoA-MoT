import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1(2) on n-dimensional complex projective space.
    """
    try:
        n_str = input("Enter the dimension 'n' of the complex projective space P^n: ")
        n = int(n_str)
        if n < 1:
            print("The dimension 'n' must be a positive integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer for 'n'.")
        return

    # h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1)) = (n+1) * (n+1)
    h0_O1_n1 = (n + 1)**2

    # h^0(P^n, O(2)) = C(n+2, 2)
    # Using integer division // for a clean result
    h0_O2 = (n + 2) * (n + 1) // 2

    # The dimension is the difference due to the short exact sequence of global sections.
    dimension = h0_O1_n1 - h0_O2

    print(f"\nFor the projective space P^{n}, the dimension of H^0(P^{n}, Omega^1(2)) is determined by the formula:")
    print(f"dim = h^0(P^{n}, O(1)^(n+1)) - h^0(P^{n}, O(2))")
    print("\nCalculating each term:")
    print(f"h^0(P^{n}, O(1)^(n+1)) = ({n}+1)^2 = {h0_O1_n1}")
    print(f"h^0(P^{n}, O(2)) = C({n}+2, 2) = ({n}+2)({n}+1)/2 = {h0_O2}")

    print("\nPlugging the numbers into the equation:")
    print(f"Dimension = {h0_O1_n1} - {h0_O2} = {dimension}")

if __name__ == '__main__':
    calculate_dimension()