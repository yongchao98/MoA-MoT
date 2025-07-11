import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_{P^n} tensored with O_{P^n}(2).
    """
    try:
        n_str = input("Please enter the value of n for the complex projective space P^n_C: ")
        n = int(n_str)
        if n < 0:
            print("The dimension n must be a non-negative integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    print(f"\nFor n = {n}:")
    print("The dimension is h^0(P^n, Omega^1(2)), which is calculated using the formula:")
    print("h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))")

    # h^0(P^n, O(1)^(n+1)) = (n+1)^2
    h0_term1 = (n + 1)**2

    # h^0(P^n, O(2)) = C(n+2, 2)
    h0_term2 = math.comb(n + 2, 2)

    # The result is the difference
    result = h0_term1 - h0_term2

    print("\nFirst, we calculate the dimensions of the related spaces:")
    print(f"h^0(P^n, O(1)^(n+1)) = (n+1)^2 = ({n}+1)^2 = {h0_term1}")
    print(f"h^0(P^n, O(2)) = C(n+2, 2) = C({n}+2, 2) = {h0_term2}")

    print("\nThen, we plug these values into the equation:")
    # Final equation with numbers
    print(f"Dimension = {h0_term1} - {h0_term2}")

    print(f"\nThe complex dimension is: {result}")

if __name__ == '__main__':
    calculate_dimension()