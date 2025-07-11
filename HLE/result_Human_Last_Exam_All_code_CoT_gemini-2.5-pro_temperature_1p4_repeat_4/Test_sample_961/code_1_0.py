# Number of groups in the free product
num_groups = 19

# Power of each commutator
power = 30

# The stable commutator length (scl) of a single commutator [a, b] in a free group F(a,b)
scl_of_base_commutator = 0.5

# Using the homogeneity property scl(g^k) = k * scl(g),
# we calculate the scl of c_i^30 in F_i.
scl_of_powered_commutator = power * scl_of_base_commutator

# Using the additivity property of scl over free products,
# the total scl is the sum of the scl of each component.
# scl(c) = sum(scl(c_i^30)) for i=1 to 19.
total_scl = num_groups * scl_of_powered_commutator

# Print the final equation and the result.
print(f"The final calculation is based on the formula: (number of groups) * (scl of one powered element).")
print(f"scl(c) = {num_groups} * {int(scl_of_powered_commutator)} = {int(total_scl)}")
