import math

def calculate_dimension(n):
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1_{P^n} twisted by O(2).

    Args:
        n (int): The dimension of the complex projective space, must be a non-negative integer.

    Returns:
        None. Prints the calculation steps and the final result.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1)) = (n+1) * (n+1)
    h0_term1 = (n + 1)**2

    # h^0(P^n, O(2)) = C(n+2, 2)
    h0_term2 = math.comb(n + 2, 2)

    # The dimension is the difference between these two values
    dimension = h0_term1 - h0_term2
    
    # The final formula is C(n+1, 2)
    final_formula_result = math.comb(n+1, 2) if n >= 1 else 0

    print(f"For the complex projective space P^{n}:")
    print("The dimension is calculated as h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))")
    print(f"h^0(P^{n}, O(1)^({n}+1)) = ({n}+1)^2 = {h0_term1}")
    print(f"h^0(P^{n}, O(2)) = C({n}+2, 2) = {h0_term2}")
    print(f"Dimension = {h0_term1} - {h0_term2} = {dimension}")
    print(f"This is equal to the formula C(n+1, 2) = {final_formula_result}")


# Example: Calculate the dimension for n=3
try:
    # You can change the value of n here to see the result for other dimensions.
    n_value = 3
    calculate_dimension(n_value)
except Exception as e:
    print(f"An error occurred: {e}")
