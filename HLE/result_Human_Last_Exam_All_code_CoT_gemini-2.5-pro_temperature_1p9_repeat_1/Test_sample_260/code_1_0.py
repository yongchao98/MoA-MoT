from fractions import Fraction

# The orders of the cyclic groups forming the free product G = pi_1(Y)
n1 = 5
n2 = 8
n3 = 2

# The number of groups in the free product
k = 3

# The index of the kernel K in G is the order of the homology group H_1(Y)
index_K = n1 * n2 * n3

# Calculate the Euler characteristic of G = pi_1(Y)
chi_G = Fraction(1, n1) + Fraction(1, n2) + Fraction(1, n3) - (k - 1)

# Calculate the Euler characteristic of the kernel K
# chi(K) = [G:K] * chi(G)
chi_K = index_K * chi_G

# The rank of a free group F_r is r = 1 - chi(F_r).
# So, the rank of K is r = 1 - chi(K).
rank_K = 1 - chi_K

# --- Output the calculation step-by-step ---
print(f"The fundamental groups are Z_{n1}, Z_{n2}, and Z_{n3}.")
print(f"The fundamental group of the connected sum is pi_1(Y) = Z_{n1} * Z_{n2} * Z_{n3}.")
print(f"The first homology group H_1(Y) is the direct sum Z_{n1} + Z_{n2} + Z_{n3}.")
print(f"The order of H_1(Y), which is the index of the kernel K, is {n1} * {n2} * {n3} = {index_K}.")
print("\nThe rank 'r' of the kernel K (a free group) is calculated using Euler characteristics:")
print(f"r = 1 - chi(K)")
print(f"r = 1 - [index(K)] * chi(pi_1(Y))")
print(f"r = 1 - {index_K} * (chi(Z_{n1}) + chi(Z_{n2}) + chi(Z_{n3}) - ({k}-1))")
print(f"r = 1 - {index_K} * (1/{n1} + 1/{n2} + 1/{n3} - {k-1})")

# Evaluate the expression inside the parenthesis
chi_G_val = 1/n1 + 1/n2 + 1/n3 - (k-1)
print(f"r = 1 - {index_K} * ({1/n1} + {1/n2} + {1/n2} - {k-1})")
print(f"r = 1 - {index_K} * ({chi_G_val})")

# Evaluate the product
chi_K_val = index_K * chi_G_val
print(f"r = 1 - ({chi_K_val})")

# Final result
rank_K_val = 1 - chi_K_val
print(f"r = {int(rank_K_val)}")

print("\nFinal Answer:")
print(f"The rank of the kernel K as a free group is {int(rank_K)}.")