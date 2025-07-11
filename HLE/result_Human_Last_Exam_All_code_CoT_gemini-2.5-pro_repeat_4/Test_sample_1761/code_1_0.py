import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_{P^n} tensored with O(2) on complex projective space P^n.

    The dimension is given by the formula: n*(n+1)/2.
    This script prompts the user for the value of n and shows the calculation
    based on the long exact sequence in cohomology derived from the Euler sequence.
    """
    try:
        n_str = input("Enter the dimension 'n' of the complex projective space P^n: ")
        n = int(n_str)
        if n < 0:
            print("Please enter a non-negative integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # Using the standard formula for cohomology of line bundles on P^n:
    # h^0(P^n, O(k)) = C(n+k, k) for k >= 0.

    # Calculate h^0(P^n, O(1))
    h0_O1 = math.comb(n + 1, 1)
    
    # Calculate h^0(P^n, O(1)^(n+1))
    h0_O1_sum = (n + 1) * h0_O1
    
    # Calculate h^0(P^n, O(2))
    h0_O2 = math.comb(n + 2, 2)

    # From the long exact sequence in cohomology, we have the relation:
    # h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))
    result = h0_O1_sum - h0_O2
    
    print(f"\nFor n = {n}, the calculation is as follows:")
    print("-" * 50)
    print("The equation for the dimension comes from the long exact sequence in cohomology:")
    print("h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))")
    print("\nCalculating each term:")
    
    # h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1))
    # h^0(P^n, O(1)) = C(n+1, 1)
    equation_part1_val = n + 1
    print(f"h^0(P^{n}, O(1)^(n+1)) = ({n}+1) * C({n}+1, 1) = {n+1} * {equation_part1_val} = {h0_O1_sum}")

    # h^0(P^n, O(2)) = C(n+2, 2)
    print(f"h^0(P^{n}, O(2)) = C({n}+2, 2) = {h0_O2}")

    print("\nSubstituting these values into the equation:")
    print(f"h^0(P^{n}, Omega^1(2)) = {h0_O1_sum} - {h0_O2} = {result}")
    
    print("-" * 50)
    print(f"The final complex dimension is {result}.")

if __name__ == "__main__":
    calculate_dimension()