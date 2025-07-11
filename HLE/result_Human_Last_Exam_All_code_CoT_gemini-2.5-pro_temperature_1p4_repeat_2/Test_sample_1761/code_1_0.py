def combinations(n, k):
    """
    Helper function to calculate combinations C(n, k) for non-negative integers.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def calculate_dimension_of_sections():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_P^n(2) on complex projective space P^n.
    """
    try:
        n_str = input("Please enter the dimension n of the complex projective space P^n (a positive integer): ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension n must be a positive integer.")
            return

        # h^0(P^n, O(1)) = C(n+1, 1)
        h0_O1 = n + 1
        
        # h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1))
        h0_O1_oplus = (n + 1) * h0_O1
        
        # h^0(P^n, O(2)) = C(n+2, 2)
        h0_O2 = combinations(n + 2, 2)
        
        # The result h^0(P^n, Omega^1(2))
        result = h0_O1_oplus - h0_O2
        
        print(f"\nTo calculate the dimension h^0(P^n, Omega^1(2)) for n = {n}, we use the formula derived from the Euler sequence:")
        print(f"h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^({n+1})) - h^0(P^n, O(2))")
        print("\nStep 1: Calculate the dimension of the constituent terms.")
        print(f"h^0(P^n, O(1)) = C({n}+1, 1) = C({n+1}, 1) = {h0_O1}")
        print(f"h^0(P^n, O(1)^({n+1})) = ({n}+1) * h^0(P^n, O(1)) = {n+1} * {h0_O1} = {h0_O1_oplus}")
        print(f"h^0(P^n, O(2)) = C({n}+2, 2) = C({n+2}, 2) = {h0_O2}")
        
        print("\nStep 2: Substitute the values into the formula.")
        print(f"Dimension = {h0_O1_oplus} - {h0_O2} = {result}")
        
        # Verification using the final simplified formula C(n+1, 2)
        direct_result = combinations(n + 1, 2)
        print(f"\nThis confirms the simplified formula: C(n+1, 2) = C({n+1}, 2) = {direct_result}.")

    except ValueError:
        print("Error: Invalid input. Please enter an integer for n.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_dimension_of_sections()