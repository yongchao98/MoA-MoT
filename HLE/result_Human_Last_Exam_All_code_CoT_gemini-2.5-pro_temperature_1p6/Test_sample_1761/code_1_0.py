import math

def h0(n, k):
    """
    Calculates the dimension of the space of global sections of O(k) on P^n,
    which is h^0(P^n, O(k)) = C(n+k, k) for k >= 0.
    """
    if k < 0:
        return 0
    # math.comb(a, b) calculates "a choose b"
    return math.comb(n + k, k)

def solve_dimension_for_n():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1(2) on P^n for a user-specified n.
    """
    try:
        # Prompt user for the dimension 'n' of the projective space
        n_str = input("Please enter the dimension n of the complex projective space P^n (a non-negative integer): ")
        n = int(n_str)
        
        if n < 0:
            print("Error: The dimension n must be a non-negative integer.")
            return

        # Calculate the dimension of the first term: h^0(O(1)^{n+1}) = (n+1)^2
        dim_term1 = (n + 1) ** 2

        # Calculate the dimension of the second term: h^0(O(2)) = C(n+2, 2)
        dim_term2 = h0(n, 2)

        # The final result is the difference
        result = dim_term1 - dim_term2
        
        # Output the calculation process, including each number in the final equation
        print(f"\nFor n = {n}, the dimension is calculated using the formula derived from the Euler sequence:")
        print(f"h^0(P^{n}, Omega^1(2)) = h^0(P^{n}, O(1)^({n+1})) - h^0(P^{n}, O(2))")
        
        # Displaying the calculation of each number
        print(f"                   = ({n}+1)^2 - C({n}+2, 2)")
        print(f"                   = {dim_term1} - {dim_term2}")
        print(f"                   = {result}")
        
        # As a check, show the result from the simplified final formula
        final_formula_result = n * (n + 1) // 2
        print(f"\nThe simplified general formula for the dimension is n(n+1)/2.")
        print(f"For n={n}, this evaluates to {n}*({n}+1)/2 = {final_formula_result}, which matches the result.")

    except ValueError:
        print("Invalid input. Please enter a non-negative integer for n.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    solve_dimension_for_n()