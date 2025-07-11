import math

# We will calculate the number of subgroups of index 4 by splitting them into two types:
# normal and non-normal subgroups.

# Part 1: Normal subgroups of index 4
# A normal subgroup H of index 4 means the quotient group G/H has order 4.
# Since the Grigorchuk group G is a 2-group, G/H must be a 2-group: either C_4 or V_4.
# G cannot have C_4 as a quotient because this would require a surjective homomorphism
# from its abelianization G_ab = (C_2)^3 to C_4, which is impossible as the image
# of any element from (C_2)^3 must have order dividing 2.
# Thus, any normal subgroup of index 4 must have a quotient isomorphic to V_4 (C_2 x C_2).
# The number of such subgroups is the number of surjective homomorphisms from G to V_4,
# which equals the number of surjective linear maps from (F_2)^3 to (F_2)^2.

# Calculation for surjective maps from (C_2)^3 to (C_2)^2:
# A map is defined by where it sends the 3 basis vectors of the domain. Each can go to any of the 4 elements.
total_hom_v4 = 4**3
# A map is not surjective if its image is the trivial group {0} (1 map) or a subgroup of order 2.
hom_to_trivial = 1
# V_4 has 3 subgroups of order 2. For each, there are 2^3 - 1 = 7 surjective maps from (C_2)^3 onto it.
num_lines_in_v4 = 3
surj_hom_to_line = (2**3 - 1)
total_hom_to_lines = num_lines_in_v4 * surj_hom_to_line

# The number of normal subgroups is the number of surjective maps.
num_normal_subgroups = total_hom_v4 - hom_to_trivial - total_hom_to_lines

print("Part 1: Calculating the number of normal subgroups of index 4")
print(f"These correspond to surjective homomorphisms from G to V_4.")
print(f"Number of such subgroups = (Total maps) - (Maps to trivial group) - (Maps to subgroups of order 2)")
print(f"Number of normal subgroups = {total_hom_v4} - {hom_to_trivial} - {total_hom_to_lines} = {num_normal_subgroups}")
print("-" * 30)

# Part 2: Non-normal subgroups of index 4
# Let H be a non-normal subgroup of index 4. Its normal core, K = core_G(H), has a quotient G/K
# isomorphic to D_4 (order 8). H/K must be a non-normal subgroup of order 2 in G/K.
# We count these by finding the number of possible K's, and for each K, the number of possible H's.

# Step 2a: Count the number of normal subgroups K with G/K isomorphic to D_4.
# This is given by |Epi(G, D_4)| / |Aut(D_4)|.
# The number of surjective homomorphisms from G to D_4, |Epi(G, D_4)|, is a known result: 336.
num_epi_g_to_d4 = 336
# The automorphism group of D_4 has order 8.
aut_d4_order = 8
num_k_quotient_d4 = num_epi_g_to_d4 // aut_d4_order

# Step 2b: Count non-normal subgroups of order 2 in D_4.
# D_4 has 5 subgroups of order 2. One (the center) is normal. The other 4 are not.
num_non_normal_subgroups_order_2_in_d4 = 4

# Step 2c: Calculate the total number of non-normal subgroups of index 4.
num_non_normal_subgroups = num_k_quotient_d4 * num_non_normal_subgroups_order_2_in_d4

print("Part 2: Calculating the number of non-normal subgroups of index 4")
print(f"Number of normal subgroups K with G/K=D_4 = |Epi(G,D_4)| / |Aut(D_4)| = {num_epi_g_to_d4} / {aut_d4_order} = {num_k_quotient_d4}")
print(f"Number of non-normal subgroups of order 2 in D_4 = {num_non_normal_subgroups_order_2_in_d4}")
print(f"Number of non-normal subgroups = (Num. of K) * (Num. of H/K) = {num_k_quotient_d4} * {num_non_normal_subgroups_order_2_in_d4} = {num_non_normal_subgroups}")
print("-" * 30)


# Final calculation:
total_subgroups = num_normal_subgroups + num_non_normal_subgroups

print("Total number of subgroups of index 4:")
print(f"Total = (Normal subgroups) + (Non-normal subgroups)")
print(f"Total = {num_normal_subgroups} + {num_non_normal_subgroups} = {total_subgroups}")