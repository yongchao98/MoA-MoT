import sympy

def calculate_dimension_formula():
    """
    Symbolically calculates the complex dimension of H^0(P^n, Omega^1(2)) for a general n.
    """
    n = sympy.Symbol('n')

    print("We want to compute the dimension h^0(Ω^1(2)) as a function of n.")
    print("The formula derived from the dual Euler sequence is:")
    print("h^0(Ω^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))\n")

    # The dimension of global sections of O(k) is the number of
    # homogeneous polynomials of degree k in n+1 variables, which is C(n+k, k).
    
    # Calculate h^0(O(1)^(n+1))
    print("Step 1: Calculate the dimension of the middle term, h^0(O(1)^(n+1)).")
    h0_O1 = sympy.binomial(n + 1, 1)
    print(f"First, h^0(O(1)) = C(n+1, 1) = {h0_O1}")
    
    h0_O1_oplus_n1 = (n + 1) * h0_O1
    # Expand the expression for better readability
    expanded_term1 = sympy.expand(h0_O1_oplus_n1)
    print(f"Then, h^0(O(1)^(n+1)) = (n+1) * h^0(O(1)) = (n+1) * ({h0_O1}) = {expanded_term1}\n")
    
    # Calculate h^0(O(2))
    print("Step 2: Calculate the dimension of the rightmost term, h^0(O(2)).")
    h0_O2 = sympy.binomial(n + 2, 2)
    # Expand the binomial coefficient
    expanded_term2 = sympy.expand(h0_O2)
    print(f"h^0(O(2)) = C(n+2, 2) = {expanded_term2}\n")
    
    # Final Calculation
    print("Step 3: Substitute these into the formula for the final equation.")
    print(f"h^0(Ω^1(2)) = {expanded_term1} - ({expanded_term2})")

    # The dimension is the difference between these two values.
    result = h0_O1_oplus_n1 - h0_O2
    
    # Simplify the final result to get a clean formula.
    simplified_result = sympy.simplify(result)
    
    print("\nSimplifying the expression gives the final formula for the dimension:")
    print(f"h^0(Ω^1(2)) = {simplified_result}")

calculate_dimension_formula()