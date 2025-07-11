# Number of free groups in the product
num_groups = 19

# Exponent of each commutator
exponent = 30

# The stable commutator length (scl) of a single commutator [a_i, b_i]
# in a free group on two generators is a known value.
scl_of_single_commutator = 0.5

# The calculation combines these values based on the properties of scl.
# scl(c) = exponent * (num_groups * scl_of_single_commutator)

# Perform the calculation
intermediate_sum = num_groups * scl_of_single_commutator
final_scl = exponent * intermediate_sum

# Print the final equation with all its numerical components, as requested.
print("The final calculation is an equation formed by the problem's parameters:")
print(f"{exponent} * ({num_groups} * {scl_of_single_commutator}) = {final_scl}")