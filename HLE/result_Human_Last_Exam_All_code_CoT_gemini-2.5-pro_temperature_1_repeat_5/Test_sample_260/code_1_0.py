import math

# Step 1: Define the orders of the fundamental groups of X1, X2, X3.
# pi_1(X1) = Z_5, pi_1(X2) = Z_8, pi_1(X3) = Z_2.
orders = [5, 8, 2]
n = len(orders)

# Step 2: Calculate the index [G:K] = |H_1(Y)|.
# H_1(Y) = Z_5 + Z_8 + Z_2.
# |H_1(Y)| = 5 * 8 * 2.
index = math.prod(orders)

# Step 3: Calculate the Euler characteristic of G = pi_1(Y).
# chi(G) = (1/|G1| + 1/|G2| + 1/|G3|) - (n-1).
sum_reciprocals = sum(1/o for o in orders)
chi_G = sum_reciprocals - (n - 1)

# Step 4: Use the formula chi(K) = [G:K] * chi(G).
chi_K = index * chi_G

# Step 5: The rank r of the free group K is given by chi(K) = 1 - r.
# So, r = 1 - chi(K).
rank = 1 - chi_K

# Print the calculation steps
print(f"The fundamental groups of the spaces are G1=Z_{orders[0]}, G2=Z_{orders[1]}, and G3=Z_{orders[2]}.")
print(f"The fundamental group of their connected sum is G = G1 * G2 * G3.")
print(f"The kernel K of the Hurewicz map is the commutator subgroup [G, G].")
print("The rank of K can be found using the formula r = 1 - [G:K] * chi(G).")
print("-" * 20)
print(f"The index of the subgroup K in G is [G:K] = |H1(Y)| = {orders[0]} * {orders[1]} * {orders[2]} = {index}")
print(f"The Euler characteristic of G is chi(G) = (1/{orders[0]} + 1/{orders[1]} + 1/{orders[2]}) - ({n}-1) = {sum_reciprocals} - {n-1} = {chi_G}")
print(f"The Euler characteristic of K is chi(K) = [G:K] * chi(G) = {index} * ({chi_G}) = {chi_K}")
print(f"The rank of K is r = 1 - chi(K) = 1 - ({chi_K}) = {int(rank)}")
print("-" * 20)
print(f"Final calculation: rank = 1 - {index} * ( (1/{orders[0]} + 1/{orders[1]} + 1/{orders[2]}) - ({n}-1) ) = {int(rank)}")
