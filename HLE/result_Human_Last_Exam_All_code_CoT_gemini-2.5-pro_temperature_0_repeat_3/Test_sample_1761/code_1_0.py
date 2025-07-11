import sympy

def solve_dimension():
    """
    Calculates the complex dimension of the space of global sections of
    the sheaf Omega^1_{P^n} tensored by O_{P^n}(2).
    """
    # Define n as a symbolic variable
    n = sympy.Symbol('n', integer=True, positive=True)

    # The dimension is given by the formula derived from the Euler sequence:
    # h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^(n+1)) - h^0(P^n, O(2))
    print("The dimension h^0(P^n, Ω^1(2)) is calculated from the exact sequence of global sections:")
    print("0 -> H^0(Ω^1(2)) -> H^0(O(1)^(n+1)) -> H^0(O(2)) -> 0")
    print("\nThis gives the relation for dimensions:")
    print("h^0(Ω^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))\n")

    # Calculate the dimension of H^0(P^n, O(1)^(n+1))
    # h^0(P^n, O(1)) = C(n+1, n)
    h0_O1 = sympy.binomial(n + 1, n)
    # h^0(P^n, O(1)^(n+1)) = (n+1) * h^0(P^n, O(1))
    h0_O1_sum = (n + 1) * h0_O1
    
    print("Step 1: Calculate h^0(O(1)^(n+1))")
    print(f"h^0(O(1)) = C(n+1, n) = {h0_O1}")
    print(f"h^0(O(1)^(n+1)) = (n+1) * h^0(O(1)) = ({n+1}) * ({h0_O1}) = {sympy.expand(h0_O1_sum)}")
    print("-" * 40)

    # Calculate the dimension of H^0(P^n, O(2))
    # h^0(P^n, O(2)) = C(n+2, n)
    h0_O2 = sympy.binomial(n + 2, n)
    
    print("Step 2: Calculate h^0(O(2))")
    print(f"h^0(O(2)) = C(n+2, n) = {h0_O2}")
    print("-" * 40)

    # Substitute these into the main formula
    print("Step 3: Substitute and simplify")
    print(f"h^0(Ω^1(2)) = ({sympy.expand(h0_O1_sum)}) - ({h0_O2})")

    # Perform the subtraction and simplify
    result = h0_O1_sum - h0_O2
    simplified_result = sympy.simplify(result)

    print("\nSimplifying the expression:")
    # To show the common factor step
    common_factor_expr = sympy.factor(result)
    print(f"= {common_factor_expr}")
    
    print("\nThe final formula for the dimension is:")
    print(simplified_result)

if __name__ == '__main__':
    solve_dimension()