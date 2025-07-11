import sympy

def calculate_dimension():
    """
    Calculates and explains the derivation of the complex dimension of the space of 
    global sections of the sheaf Omega^1_{P^n_C} tensored by O_{P^n_C}(2).
    """
    n = sympy.symbols('n')

    print("Step 1: Start with the Euler sequence on the complex projective space P^n.")
    print("0 -> Omega^1 -> O(-1)^(n+1) -> O -> 0\n")

    print("Step 2: Twist the sequence by O(2) to get a short exact sequence involving Omega^1(2).")
    print("0 -> Omega^1(2) -> O(1)^(n+1) -> O(2) -> 0\n")

    print("Step 3: From the long exact sequence in cohomology, we derive the relation:")
    print("h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))\n")

    print("Step 4: Calculate the dimensions of the terms on the right-hand side.")
    print("The dimension of global sections of O(k) on P^n is h^0(P^n, O(k)) = C(n+k, k).\n")

    # Calculate h^0(O(1)) and then h^0(O(1)^(n+1))
    h0_O1 = sympy.binomial(n + 1, 1)
    h0_O1_oplus_n1 = (n + 1) * h0_O1
    
    # Calculate h^0(O(2))
    h0_O2 = sympy.binomial(n + 2, 2)

    print(f"The dimension of H^0(P^n, O(1)^(n+1)) is:")
    print(f"h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1)) = (n+1) * C(n+1, 1) = (n+1) * {h0_O1} = {sympy.expand(h0_O1_oplus_n1)}\n")
    
    print(f"The dimension of H^0(P^n, O(2)) is:")
    print(f"h^0(P^n, O(2)) = C(n+2, 2) = {h0_O2}\n")

    print("Step 5: Substitute these into the equation to find the final dimension.")
    # The final equation with each 'number' (symbolic expression) shown
    print("The equation is:")
    print(f"h^0(P^n, Omega^1(2)) = {sympy.expand(h0_O1_oplus_n1)} - ({h0_O2})")

    # Simplify the final expression
    final_dimension = sympy.simplify(h0_O1_oplus_n1 - h0_O2)

    print("\nSimplifying the expression leads to the final result:")
    print(f"h^0(P^n, Omega^1(2)) = {final_dimension}\n")
    
    print("This can also be written as a binomial coefficient:")
    print(f"h^0(P^n, Omega^1(2)) = C(n+1, 2)")


if __name__ == '__main__':
    calculate_dimension()