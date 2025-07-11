import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_P^n(2) on complex projective space P^n.

    The problem is solved for the case n=2.
    """
    # The dimension 'n' of the projective space. We choose n=2 as an example.
    n = 2

    # The dimension is given by the formula:
    # h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))

    # Calculate h^0(O(1)^(n+1)) = (n+1) * C(n+1, 1) = (n+1)^2
    h0_O1_term = (n + 1)**2

    # Calculate h^0(O(2)) = C(n+2, 2)
    h0_O2_term = math.comb(n + 2, 2)

    # The final dimension is the difference.
    dimension = h0_O1_term - h0_O2_term
    
    print(f"This script calculates the dimension for n = {n}.")
    print("The formula for the dimension is h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^({n+1})) - h^0(P^n, O(2)).")
    print(f"For n={n}, the first term is h^0(P^{n}, O(1)^({n+1})) = {h0_O1_term}.")
    print(f"For n={n}, the second term is h^0(P^{n}, O(2)) = {h0_O2_term}.")
    print(f"The final equation with these numbers is: {h0_O1_term} - {h0_O2_term} = {dimension}")
    print(f"\nThe complex dimension is {dimension}.")

calculate_dimension()