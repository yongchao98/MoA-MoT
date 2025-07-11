import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1_P^n tensored with O(2) on complex projective space P^n.

    The user is prompted to enter the dimension 'n' of the projective space.
    The script then prints the steps of the calculation based on the derived formula.
    """
    try:
        n_str = input("Enter the dimension 'n' of the complex projective space P^n: ")
        n = int(n_str)
        if n < 1:
            print("Error: The dimension 'n' must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # From the long exact sequence of cohomology, we have the formula:
    # h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))

    # Calculate h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1)) = (n+1) * (n+1)
    h0_term1 = (n + 1)**2

    # Calculate h^0(P^n, O(2)) = C(n+2, 2)
    h0_term2 = math.comb(n + 2, 2)

    # Calculate the final dimension
    result = h0_term1 - h0_term2

    print("\n" + "="*50)
    print(f"Calculating the dimension of H^0(P^{n}, Omega^1(2)) for n = {n}")
    print("="*50)
    print("\nThe dimension is computed using the formula:")
    print("h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))")
    
    print("\nFirst, we calculate the dimension of each term:")
    print(f"h^0(P^{n}, O(1)^(n+1)) = (n+1)^2 = ({n}+1)^2 = {h0_term1}")
    print(f"h^0(P^{n}, O(2)) = C(n+2, 2) = C({n}+2, 2) = {h0_term2}")

    print("\nNow, we substitute these values into the formula:")
    print(f"Dimension = {h0_term1} - {h0_term2}")
    
    print(f"\nThe complex dimension of the space of global sections is: {result}")
    
    # Verify with the simplified formula C(n+1, 2)
    final_formula_result = math.comb(n + 1, 2)
    print(f"\nNote: This result is equivalent to the binomial coefficient C(n+1, 2) = C({n}+1, 2) = {final_formula_result}.")

if __name__ == '__main__':
    calculate_dimension()