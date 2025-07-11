import math
from fractions import Fraction

# This script calculates the rank of the kernel of the Hurewicz homomorphism
# for the connected sum of three specific topological spaces.

# --- Explanation of the Method ---
print("This script solves the problem by following these steps from algebraic topology:")
print("1. The fundamental group of the connected sum Y is the free product G = Z_5 * Z_8 * Z_2.")
print("2. The kernel K of the Hurewicz map is the commutator subgroup [G, G], which is a free group.")
print("3. The rank 'r' of K is found using the Euler characteristic formula: r = 1 - chi(K).")
print("4. The Euler characteristic of K is computed via the relation chi(K) = [G:K] * chi(G), where [G:K] is the index of K in G.")
print("\n--- Calculation Steps ---")

# --- Step 1: Define group properties ---
# The orders of the cyclic groups from the fundamental groups of X1, X2, and X3.
orders = [5, 8, 2]
num_groups = len(orders)
print(f"The fundamental group is G = Z_{orders[0]} * Z_{orders[1]} * Z_{orders[2]}.")

# --- Step 2: Calculate the Euler characteristic of G ---
# chi(G) = sum(chi(Z_n)) - (k-1) = sum(1/n) - (k-1)
chi_G_terms_str = [f"1/{n}" for n in orders]
chi_G = sum(Fraction(1, n) for n in orders) - (num_groups - 1)
print(f"\nThe Euler characteristic of G is calculated as:")
print(f"chi(G) = ({' + '.join(chi_G_terms_str)}) - ({num_groups} - 1) = {chi_G.numerator}/{chi_G.denominator}")

# --- Step 3: Calculate the index of the kernel K in G ---
# The index [G:K] is the order of the abelianization G/[G,G] = H_1(Y).
# |H_1(Y)| = |Z_5 + Z_8 + Z_2| = 5 * 8 * 2.
index = math.prod(orders)
index_calc_str = ' * '.join(map(str, orders))
print(f"\nThe index [G:K] is the order of the first homology group H_1(Y).")
print(f"[G:K] = |H_1(Y)| = {index_calc_str} = {index}")

# --- Step 4: Calculate the Euler characteristic of K ---
# chi(K) = [G:K] * chi(G)
chi_K = index * chi_G
print(f"\nThe Euler characteristic of the kernel K is:")
print(f"chi(K) = [G:K] * chi(G) = {index} * ({chi_G.numerator}/{chi_G.denominator}) = {int(chi_K)}")

# --- Step 5: Calculate the rank of K ---
# For a free group K of rank r, chi(K) = 1 - r.
# So, r = 1 - chi(K).
rank = 1 - chi_K
print(f"\nThe rank 'r' of the free group K is found using the formula r = 1 - chi(K).")
print("The final equation for the rank is:")
# The problem asks to output each number in the final equation.
# The numbers are 1 and chi_K.
print(f"r = 1 - ({int(chi_K)}) = {int(rank)}")
