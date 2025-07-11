def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1(2) on n-dimensional complex projective space.
    """
    try:
        n_str = input("Enter the dimension 'n' of the complex projective space P^n (n >= 1): ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension 'n' must be a positive integer.")
            return
    except ValueError:
        print("Error: Please enter a valid integer for 'n'.")
        return

    # h^0(P^n, O(1)^(n+1)) = (n+1)^2
    h0_O1_oplus_n1 = (n + 1)**2

    # h^0(P^n, O(2)) = C(n+2, 2) = (n+2)(n+1)/2
    h0_O2 = (n + 2) * (n + 1) // 2
    
    # The final dimension
    dimension = h0_O1_oplus_n1 - h0_O2
    
    print("\n" + "="*40)
    print(f"Calculating dimension for P^{n}:")
    print("The dimension is h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))")
    print(f"h^0(P^{n}, O(1)^({n+1})) = ({n}+1)^2 = {h0_O1_oplus_n1}")
    print(f"h^0(P^{n}, O(2)) = C({n}+2, 2) = {h0_O2}")
    print(f"\nFinal Dimension = {h0_O1_oplus_n1} - {h0_O2} = {dimension}")
    print("="*40)

if __name__ == "__main__":
    calculate_dimension()