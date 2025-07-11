# Define the numbers from the problem description
num_groups = 19
exponent = 30
scl_of_base_commutator = 0.5

# The calculation is based on the properties of stable commutator length (scl).
# scl(c) = scl(product_{i=1 to 19} c_i^30)
# By additivity over free products:
# scl(c) = sum_{i=1 to 19} scl(c_i^30)
# By homogeneity:
# scl(c) = sum_{i=1 to 19} 30 * scl(c_i)
# Given scl(c_i) = 1/2 for c_i = [a_i, b_i] in a free group F_2:
# scl(c) = 19 * 30 * (1/2)

# Perform the final calculation
total_scl = num_groups * exponent * scl_of_base_commutator

# Output the equation with the numbers and the final result
print(f"The final calculation is:")
print(f"{num_groups} * {exponent} * {scl_of_base_commutator} = {total_scl}")