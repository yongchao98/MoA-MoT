# Number of free groups in the free product
num_groups = 19

# The exponent for each commutator
exponent = 30

# The stable commutator length of a single commutator c_i = [a_i, b_i]
# in the free group F_i = <a_i, b_i> is a known value of 1/2.
scl_of_one_commutator = 0.5

# The total stable commutator length is calculated by applying the properties of scl:
# 1. Additivity over free products: scl(product) = sum(scl of parts)
# 2. Homogeneity: scl(g^k) = k * scl(g)
#
# scl(c) = scl(c_1^30 * ... * c_19^30)
#        = sum_{i=1 to 19} scl(c_i^30)
#        = sum_{i=1 to 19} 30 * scl(c_i)
#        = sum_{i=1 to 19} 30 * 0.5
#        = 19 * 30 * 0.5

total_scl = int(num_groups * exponent * scl_of_one_commutator)

# Print the final equation showing each number involved in the calculation
print(f"The stable commutator length is calculated as:")
print(f"{num_groups} * {exponent} * {scl_of_one_commutator} = {total_scl}")