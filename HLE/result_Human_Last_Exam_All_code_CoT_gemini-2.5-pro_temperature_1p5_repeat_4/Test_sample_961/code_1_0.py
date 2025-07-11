# The number of free groups in the free product.
num_groups = 19

# The exponent for each commutator c_i.
exponent = 30

# The stable commutator length (scl) of a single commutator [a, b] in a free group F(a,b)
# is a known theoretical result, equal to 1/2.
scl_of_basic_commutator = 0.5

# Step 1: According to the homogeneity property of scl, scl(g^n) = n * scl(g).
# We calculate the scl for a single component c_i^30 in its respective group F_i.
# scl_{F_i}(c_i^30) = 30 * scl_{F_i}(c_i) = 30 * 0.5 = 15.
scl_of_one_powered_component = exponent * scl_of_basic_commutator

# Step 2: Due to the additivity of scl over free products, the total scl is the sum
# of the scls of the components. As all components are structurally identical,
# this simplifies to a multiplication.
# scl_G(c) = sum_{i=1 to 19} scl_{F_i}(c_i^30) = 19 * 15.
total_scl = num_groups * scl_of_one_powered_component

# Print the final equation with all its numbers, as requested.
print(f"The final calculation for the stable commutator length is:")
print(f"{num_groups} * {int(scl_of_one_powered_component)} = {int(total_scl)}")