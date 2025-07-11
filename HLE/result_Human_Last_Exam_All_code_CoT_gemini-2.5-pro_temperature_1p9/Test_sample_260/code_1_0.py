from fractions import Fraction

# Define the orders of the cyclic groups from the problem description.
# pi_1(X_1) is Z_5, pi_1(X_2) is Z_8, pi_1(X_3) is Z_2.
n1, n2, n3 = 5, 8, 2
orders = [n1, n2, n3]
num_groups = len(orders)

# Let G = pi_1(Y) be the free product of these groups. The kernel K of the
# Hurewicz map G -> H_1(Y) is the commutator subgroup [G, G].
# K is a free group. Its rank 'r' is determined by its Euler characteristic chi(K) = 1 - r.
# We use the formula: chi(K) = [G:K] * chi(G).

# Step 1: Calculate the index [G:K]. This is the order of the first homology group H_1(Y).
# H_1(Y) = Z_5 x Z_8 x Z_2
index_K = n1 * n2 * n3

# Step 2: Calculate the Euler characteristic of G.
# chi(G) = sum(chi(Z_n)) - (number of groups - 1) where chi(Z_n) = 1/n.
sum_chi_components = sum(Fraction(1, n) for n in orders)
chi_G = sum_chi_components - (num_groups - 1)

# Step 3: Calculate the Euler characteristic of the kernel K.
chi_K = index_K * chi_G
# Since K is a free group, chi(K) must be an integer.
chi_K_int = int(chi_K)

# Step 4: Calculate the rank r from chi(K) = 1 - r.
rank_K = 1 - chi_K_int

# Display the step-by-step calculation with the specific numbers.
print(f"The rank 'r' of the kernel K is calculated as r = 1 - chi(K).")
print(f"First, we find chi(K) using the formula chi(K) = [G:K] * chi(G).\n")

print(f"1. The index [G:K] is the order of H_1(Y) = Z_{n1} x Z_{n2} x Z_{n3}:")
print(f"   [G:K] = {n1} * {n2} * {n3} = {index_K}.\n")

print(f"2. The Euler characteristic chi(G) is calculated from its components:")
print(f"   chi(G) = (1/{n1} + 1/{n2} + 1/{n3}) - ({num_groups} - 1) = {chi_G}.\n")

print(f"3. Using these, the Euler characteristic of the kernel is:")
print(f"   chi(K) = {index_K} * ({chi_G}) = {chi_K_int}.\n")

print(f"4. Finally, the rank is derived from chi(K):")
print(f"   r = 1 - chi(K) = 1 - ({chi_K_int}) = {rank_K}.\n")

print("The complete equation for the rank is:")
print(f"r = 1 - ({index_K}) * ( (1/{n1} + 1/{n2} + 1/{n3}) - ({num_groups - 1}) ) = {rank_K}")