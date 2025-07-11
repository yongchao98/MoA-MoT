from fractions import Fraction

# Step 1: Define the orders of the fundamental groups of X_1, X_2, and X_3.
# |pi_1(X_1)| = 5
# |pi_1(X_2)| = 8
# |pi_1(X_3)| = 2
n1, n2, n3 = 5, 8, 2
orders = [n1, n2, n3]
k = len(orders)

print(f"The orders of the fundamental groups are {n1}, {n2}, and {n3}.")

# Step 2: Calculate the index of the kernel K in G = pi_1(Y).
# This is the order of the first homology group H_1(Y), which is the abelianization of G.
# [G:K] = |H_1(Y)| = |Z_5 x Z_8 x Z_2| = 5 * 8 * 2
index = n1 * n2 * n3
print(f"The index [G:K] is the product of the orders: {n1} * {n2} * {n3} = {index}")

# Step 3: Calculate the Euler characteristic of the group G = pi_1(Y).
# chi(G) = (1/|G_1| + 1/|G_2| + ... + 1/|G_k|) - (k-1)
sum_recip = sum(Fraction(1, n) for n in orders)
chi_G = sum_recip - (k - 1)
print(f"The Euler characteristic of G is chi(G) = (1/{n1} + 1/{n2} + 1/{n3}) - ({k}-1) = {sum_recip} - {k-1} = {chi_G}")

# Step 4: Calculate the Euler characteristic of the kernel K using the formula chi(K) = [G:K] * chi(G).
# K is a free group of rank r, so chi(K) = 1 - r.
chi_K = index * chi_G
print(f"The Euler characteristic of the kernel K is chi(K) = {index} * ({chi_G}) = {int(chi_K)}")

# Step 5: Solve for the rank r from the equation 1 - r = chi(K).
r = 1 - chi_K
print(f"\nThe rank r of the kernel K is found by solving 1 - r = chi(K).")
print(f"The final equation for the rank r is:")
print(f"r = 1 - ({int(chi_K)}) = {int(r)}")