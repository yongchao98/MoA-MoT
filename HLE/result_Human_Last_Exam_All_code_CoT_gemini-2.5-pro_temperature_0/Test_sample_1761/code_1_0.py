import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_{P^n} tensored by O_{P^n}(2).
    """
    try:
        n_str = input("Enter the dimension n of the projective space P^n (a positive integer): ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension n must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # The dimension is given by the formula: h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))
    
    # Calculate h^0(O(1)^(n+1)) = (n+1) * h^0(O(1))
    # h^0(O(1)) = C(n+1, 1) = n+1
    h0_O1 = n + 1
    h0_O1_oplus_n_plus_1 = (n + 1) * h0_O1

    # Calculate h^0(O(2)) = C(n+2, 2)
    h0_O2 = math.comb(n + 2, 2)

    # The final dimension is the difference
    dimension = h0_O1_oplus_n_plus_1 - h0_O2

    print("\n--- Calculation Steps ---")
    print(f"The dimension is calculated using the formula derived from the Euler sequence:")
    print(f"h^0(P^{n}, Omega^1(2)) = h^0(P^{n}, O(1)^({n}+1)) - h^0(P^{n}, O(2))")
    
    print("\n1. Calculate h^0(P^{n}, O(1)^({n}+1)):")
    print(f"   h^0(P^{n}, O(1)) = C({n}+1, 1) = {h0_O1}")
    print(f"   h^0(P^{n}, O(1)^({n}+1)) = ({n}+1) * {h0_O1} = {h0_O1_oplus_n_plus_1}")

    print("\n2. Calculate h^0(P^{n}, O(2)):")
    print(f"   h^0(P^{n}, O(2)) = C({n}+2, 2) = {h0_O2}")

    print("\n3. Final Dimension:")
    print(f"   h^0(P^{n}, Omega^1(2)) = {h0_O1_oplus_n_plus_1} - {h0_O2} = {dimension}")
    
    # Verify with the simplified formula n*(n+1)/2
    formula_result = n * (n + 1) // 2
    print(f"\nThe final formula for the dimension is n*(n+1)/2.")
    print(f"For n={n}, this is {n}*({n}+1)/2 = {formula_result}.")


if __name__ == '__main__':
    calculate_dimension()