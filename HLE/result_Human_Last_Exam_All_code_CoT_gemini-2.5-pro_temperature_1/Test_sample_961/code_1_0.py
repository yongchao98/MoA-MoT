# Number of free groups in the free product
num_groups = 19

# The power to which each commutator is raised
power = 30

# The stable commutator length of a single commutator [a, b] in a free group is 1/2.
scl_of_single_commutator = 0.5

# First, we calculate the scl for one of the elements, c_i^30.
# Using the homogeneity property: scl(g^n) = n * scl(g).
scl_of_powered_commutator = power * scl_of_single_commutator

# Next, we calculate the total scl.
# Using the additivity property for free products: scl(g1*g2*...) = scl(g1) + scl(g2) + ...
total_scl = num_groups * scl_of_powered_commutator

# The final equation is the sum of the scl of each component, which simplifies to
# the number of groups multiplied by the scl of one component.
# The scl of one component c_i^30 is 30 * (1/2) = 15.
# The total scl is 19 * 15.
print(f"The final calculation is:")
print(f"{num_groups} * ({power} * {scl_of_single_commutator}) = {num_groups} * {int(scl_of_powered_commutator)} = {int(total_scl)}")
