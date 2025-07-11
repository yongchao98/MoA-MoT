# This script calculates the number of subgroups of index 4 in the Grigorchuk group.
# The calculation is based on group-theoretic properties of this group.

# Case 1: Subgroups corresponding to an action with image V_4 (Klein four-group).
# This is equivalent to the number of normal subgroups of index 4.
# This equals the number of subgroups of order 2 in the abelianization (Z/2Z)^3.
num_elements = 2**3
# Every non-identity element has order 2 and generates a unique subgroup of order 2.
num_normal_subgroups = num_elements - 1
print(f"Number of normal subgroups of index 4 (V_4 type): {num_normal_subgroups}")

# Case 2: Subgroups corresponding to an action with image D_4 (dihedral group of order 8).
# This equals the number of non-normal subgroups of order 2 in D_4, as there is a unique D_4 quotient of the Grigorchuk group.
# D_4 has 5 elements of order 2. One generates the center (a normal subgroup), and 4 generate non-normal subgroups.
num_non_normal_subgroups = 4
print(f"Number of non-normal subgroups of index 4 (D_4 type): {num_non_normal_subgroups}")

# The total number is the sum of the two disjoint cases.
total = num_normal_subgroups + num_non_normal_subgroups

print("\nThe total number of subgroups of index 4 is the sum of these two cases:")
print(f"{num_normal_subgroups} + {num_non_normal_subgroups} = {total}")