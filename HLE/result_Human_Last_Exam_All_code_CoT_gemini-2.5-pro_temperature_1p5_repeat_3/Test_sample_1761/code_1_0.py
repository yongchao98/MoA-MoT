import sympy

def solve_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_P^n(2) on complex projective n-space.
    """
    # Define n as a symbol for symbolic computation
    n = sympy.Symbol('n', integer=True, positive=True)

    # The dimension of the space of global sections H^0(P^n, O(k)) is
    # the number of monomials of degree k in n+1 variables, which is
    # given by the binomial coefficient C(n+k, k).
    def h0_Ok(n_sym, k):
        return sympy.binomial(n_sym + k, k)

    # --- Calculation Steps ---

    # We need to find h^0(Omega^1(2)). From the Euler sequence, we have the relation:
    # h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))

    # 1. Calculate h^0(O(1)^(n+1))
    # h^0(O(1)) is the dimension of the space of linear homogeneous polynomials.
    dim_H0_O1 = h0_Ok(n, 1)
    # The space O(1)^(n+1) is a direct sum of n+1 copies of O(1).
    dim_H0_O1_bundle = (n + 1) * dim_H0_O1

    # 2. Calculate h^0(O(2))
    # h^0(O(2)) is the dimension of the space of quadratic homogeneous polynomials.
    dim_H0_O2 = h0_Ok(n, 2)

    # 3. The final dimension is the difference.
    result_dim = dim_H0_O1_bundle - dim_H0_O2
    
    # --- Output the Explanation and Result ---
    
    print("This script calculates the complex dimension h^0(P^n, Omega^1(2)).")
    print("The calculation is based on the Euler sequence and sheaf cohomology.")
    print("The main formula is: h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))\n")

    print("Step 1: Calculating h^0(O(1)^(n+1))")
    print(f"h^0(O(1)) = binomial(n+1, 1) = {dim_H0_O1}")
    print(f"h^0(O(1)^(n+1)) = (n+1) * h^0(O(1)) = ({n+1}) * ({dim_H0_O1})")
    
    # Let's expand dim_H0_O1_bundle for clearer output
    expanded_dim_H0_O1_bundle = sympy.expand(dim_H0_O1_bundle)
    print(f"               = {expanded_dim_H0_O1_bundle}\n")
    
    print("Step 2: Calculating h^0(O(2))")
    print(f"h^0(O(2)) = binomial(n+2, 2) = {dim_H0_O2}\n")
    
    print("Step 3: Calculating the final dimension by subtraction.")
    print("h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))")
    # This fulfills the instruction: "output each number in the final equation!"
    print(f"               = ({expanded_dim_H0_O1_bundle}) - ({dim_H0_O2})")

    # Simplify the final expression
    simplified_result = sympy.simplify(result_dim)

    print("\nSimplifying the expression, the final dimension is:")
    print(f"h^0(Omega^1(2)) = {simplified_result}")

if __name__ == '__main__':
    solve_dimension()