import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1_P^n(C) tensored by O_P^n(C)(2) for a given n.
    """
    try:
        n_str = input("Please enter the dimension n of the complex projective space P^n_C: ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension n must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer for n.")
        return

    # Calculate h^0(P^n, O(1)^(n+1))
    # h^0(P^n, O(1)) = n+1
    h0_term1 = (n + 1)**2

    # Calculate h^0(P^n, O(2))
    # h^0(P^n, O(2)) = C(n+2, 2)
    h0_term2 = math.comb(n + 2, 2)
    
    # The dimension is the difference
    dimension = h0_term1 - h0_term2

    print(f"\nFor n = {n}, the dimension is h^0(P^{n}, Omega^1(2)).")
    print("This is computed from the exact sequence of global sections as:")
    print(f"h^0(P^{n}, Omega^1(2)) = h^0(P^{n}, O(1)^({n}+1)) - h^0(P^{n}, O(2))")
    
    print("\nCalculating the terms:")
    print(f"h^0(P^{n}, O(1)^({n}+1)) = ({n}+1)^2 = {h0_term1}")
    print(f"h^0(P^{n}, O(2)) = C({n}+2, 2) = {h0_term2}")
    
    print("\nThe final equation with the computed values is:")
    print(f"{h0_term1} - {h0_term2} = {dimension}")

if __name__ == '__main__':
    calculate_dimension()