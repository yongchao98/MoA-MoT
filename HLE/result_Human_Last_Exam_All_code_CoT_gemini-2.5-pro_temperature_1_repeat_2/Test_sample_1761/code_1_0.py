import sympy

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_P^n tensor O_P^n(2) on complex projective n-space.
    """
    # Define n as a symbolic variable for the dimension of the projective space.
    n = sympy.Symbol('n', integer=True, positive=True)

    print("This script calculates h^0(P^n, Omega^1(2)), the dimension of the space of global sections of Omega^1(2).")
    print("-" * 80)
    print("Step 1: Start with the Euler sequence for the cotangent bundle Omega^1 on P^n.")
    print("The sequence is: 0 -> Omega^1 -> O(-1)^(n+1) -> O -> 0")
    print("\nStep 2: Twist the sequence by O(2) to get a new short exact sequence.")
    print("The twisted sequence is: 0 -> Omega^1(2) -> O(1)^(n+1) -> O(2) -> 0")
    print("\nStep 3: This gives rise to a long exact sequence in cohomology:")
    print("0 -> H^0(P^n, Omega^1(2)) -> H^0(P^n, O(1)^(n+1)) -> H^0(P^n, O(2)) -> ...")
    print("-" * 80)

    print("Step 4: Calculate the dimensions of the known terms in the sequence.")
    print("We use the formula: h^0(P^n, O(k)) = C(n+k, k)")
    
    # Calculate h^0(P^n, O(1)^(n+1))
    h0_O1 = sympy.binomial(n + 1, 1)
    h0_term1 = (n + 1) * h0_O1
    
    # Calculate h^0(P^n, O(2))
    h0_term2 = sympy.binomial(n + 2, 2)
    
    print("\nDimension of the middle term H^0(P^n, O(1)^(n+1)):")
    print(f"h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1)) = (n+1) * C(n+1, 1)")
    print(f"                       = (n+1) * ({h0_O1}) = {sympy.simplify(h0_term1)}")
    
    print("\nDimension of the right term H^0(P^n, O(2)):")
    print(f"h^0(P^n, O(2)) = C(n+2, 2) = {sympy.simplify(h0_term2)}")
    print("-" * 80)

    print("Step 5: Determine the final dimension.")
    print("The map H^0(O(1)^(n+1)) -> H^0(O(2)) is surjective. This means the long exact sequence gives a short exact sequence of vector spaces:")
    print("0 -> H^0(Omega^1(2)) -> H^0(O(1)^(n+1)) -> H^0(O(2)) -> 0")
    print("\nBy the property of short exact sequences, the dimension of the first term is the dimension of the second minus the dimension of the third.")
    
    # Final calculation
    final_dimension = sympy.simplify(h0_term1 - h0_term2)
    
    print("\nThe final equation for the dimension h^0(P^n, Omega^1(2)) is:")
    print(f"h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))")
    # Outputting each "number" (symbolic part) of the final equation
    print(f"               = {sympy.simplify(h0_term1)} - ({sympy.simplify(h0_term2)})")
    print(f"               = (n+1)^2 - (n+1)*(n+2)/2")
    print(f"               = (n+1) * [ (n+1) - (n+2)/2 ]")
    print(f"               = (n+1) * [ (2n+2 - n-2)/2 ]")
    print(f"               = (n+1) * n / 2")

    print("-" * 80)
    print(f"The complex dimension is: {final_dimension}")

if __name__ == '__main__':
    calculate_dimension()