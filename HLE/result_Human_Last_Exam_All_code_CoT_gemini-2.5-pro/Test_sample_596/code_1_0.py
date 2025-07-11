# This script calculates the number of higher dimensional rooted forests on the standard
# triangulation of the Möbius band that fail to simplicially collapse onto their root set.

# Step 1: Define the dimension of the simplicial complex.
# The Möbius band is a 2D surface, so its triangulation is a 2-dimensional complex.
d = 2

# Step 2: Determine the orders of the torsion subgroups of the relevant homology groups.
# The number of non-collapsing rooted forests is given by the formula:
# N = |Tors(H_{d-1}(M; Z))| * |Tors(H_{d-2}(M; Z))| - 1
# where M is the Möbius band and Z are the integers.

# For d=2, we need H_1(M; Z) and H_0(M; Z).

# The homology groups of the Möbius band M are well-known:
# H_1(M; Z) is the group of integers, Z.
# H_0(M; Z) is the group of integers, Z (since it's path-connected).

# The torsion subgroup of Z is the trivial group {0}, which has order 1.
order_tors_h1 = 1
order_tors_h0 = 1

# Step 3: Calculate the total number of non-collapsing forests.
# The order of the group of non-collapsing forests (plus the trivial one) is the product of the orders.
group_order = order_tors_h1 * order_tors_h0

# The number of non-collapsing forests is the order of this group minus 1 (for the trivial empty forest).
num_non_collapsing_forests = group_order - 1

# Step 4: Print the final equation and the result.
print("The calculation for the number of non-collapsing rooted forests is based on the formula:")
print("N = |Tors(H_1(M))| * |Tors(H_0(M))| - 1")
print(f"The order of the torsion subgroup of H_1(M) is {order_tors_h1}.")
print(f"The order of the torsion subgroup of H_0(M) is {order_tors_h0}.")
print("\nFinal calculation:")
print(f"{order_tors_h1} * {order_tors_h0} - 1 = {num_non_collapsing_forests}")

print(f"\nThe number of such forests is: {num_non_collapsing_forests}")
