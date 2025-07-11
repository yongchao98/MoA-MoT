import math

def h0_O(n, k):
    """Calculates h^0(P^n, O(k))"""
    if k < 0:
        return 0
    # math.comb(n+k, k)
    return math.comb(n + k, k)

def solve_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_P^n tensored with O_P^n(2).
    """
    try:
        n_str = input("Please enter the dimension 'n' of the projective space P^n: ")
        n = int(n_str)
        if n < 0:
            print("Dimension 'n' must be a non-negative integer.")
            return

        # h^0(P^n, O(1)^{n+1}) = (n+1) * h^0(P^n, O(1))
        h0_O1_n_plus_1 = (n + 1) * h0_O(n, 1)

        # h^0(P^n, O(2))
        h0_O2 = h0_O(n, 2)
        
        # The result from the long exact sequence
        result = h0_O1_n_plus_1 - h0_O2
        
        # We can also compute it directly from the derived formula n(n+1)/2
        formula_result = n * (n + 1) // 2

        print("\nStep-by-step calculation:")
        print(f"h^0(P^{n}, O(1)) = {h0_O(n, 1)}")
        print(f"h^0(P^{n}, O(1)^(n+1)) = (n+1) * h^0(P^{n}, O(1)) = {n+1} * {h0_O(n, 1)} = {h0_O1_n_plus_1}")
        print(f"h^0(P^{n}, O(2)) = {h0_O2}")
        print("\nFrom the exact sequence, the dimension is h^0(O(1)^(n+1)) - h^0(O(2)):")
        
        print(f"\nFinal Equation:")
        print(f"{h0_O1_n_plus_1} - {h0_O2} = {result}")

        print(f"\nThe complex dimension is {result}.")
        
        # Verification with the final formula
        # print(f"\nVerification with formula n(n+1)/2: {formula_result}")
        # assert result == formula_result

    except ValueError:
        print("Invalid input. Please enter an integer.")

if __name__ == '__main__':
    solve_dimension()