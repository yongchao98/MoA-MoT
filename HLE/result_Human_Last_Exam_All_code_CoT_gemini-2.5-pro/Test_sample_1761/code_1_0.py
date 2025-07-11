import math

def calculate_dimension(n):
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_P^n(2) on complex projective n-space.

    The dimension is given by the formula: h^0(Omega^1(2)) = C(n+1, 2).
    This is derived from the Euler sequence, leading to the calculation:
    h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))
                   = (n+1)^2 - C(n+2, 2)
    """
    if n < 1:
        print("The dimension n of projective space must be a positive integer.")
        return

    # h^0(O(1)^(n+1)) = (n+1) * h^0(O(1)) = (n+1) * C(n+1, 1) = (n+1)^2
    h0_term1 = (n + 1)**2

    # h^0(O(2)) = C(n+2, 2)
    h0_term2 = math.comb(n + 2, 2)

    # The result is the difference
    result = h0_term1 - h0_term2
    
    # Final check using the simplified formula C(n+1, 2)
    final_result_formula = math.comb(n + 1, 2) if n >= 1 else 0
    assert result == final_result_formula

    print(f"For the complex projective space P^n with n = {n}:")
    print("\nThe dimension is calculated as h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))")
    print("\nStep 1: Calculate h^0(O(1)^(n+1))")
    print(f"h^0(O(1)^(n+1)) = ({n}+1)^2 = {h0_term1}")
    
    print("\nStep 2: Calculate h^0(O(2))")
    print(f"h^0(O(2)) = C({n}+2, 2) = {h0_term2}")
    
    print("\nStep 3: Calculate the final dimension")
    print(f"h^0(Omega^1(2)) = {h0_term1} - {h0_term2} = {result}")

    print(f"\nThis is equal to the simplified formula C(n+1, 2) = C({n}+1, 2) = {final_result_formula}.")


# Example calculation for n=3
calculate_dimension(3)