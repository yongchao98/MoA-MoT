# Step 1: Define the constants from the problem description.

# Number of free groups in the free product, which is 19.
num_groups = 19

# The exponent of each commutator, which is 30.
exponent = 30

# The stable commutator length of a single commutator [a_i, b_i] in a free group F_i.
# This is a standard result in geometric group theory.
scl_of_one_commutator = 0.5

# Step 2: Calculate the total stable commutator length.
# We use the additivity property of scl over free products and the homogeneity property.
# scl(c) = scl(c_1^30 * ... * c_19^30) = sum_{i=1 to 19} scl(c_i^30)
#          = sum_{i=1 to 19} 30 * scl(c_i)
#          = 19 * 30 * scl([a,b])
#          = 19 * 30 * (1/2)

total_scl = num_groups * exponent * scl_of_one_commutator

# Step 3: Print the final equation and the result.
# The user wants to see the numbers in the final equation.
# We cast the final result to an integer as it is a whole number.
print("The stable commutator length of c is computed as:")
print(f"scl(c) = {num_groups} * {exponent} * {scl_of_one_commutator}")
print(f"scl(c) = {int(total_scl)}")
