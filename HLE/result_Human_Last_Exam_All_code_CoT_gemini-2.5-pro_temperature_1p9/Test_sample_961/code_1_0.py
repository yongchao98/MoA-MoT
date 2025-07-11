# Number of groups in the free product
num_groups = 19

# Exponent of each commutator
exponent = 30

# scl of a single commutator [a_i, b_i] in the free group F_i
scl_base = 0.5

# Calculate the total stable commutator length
# scl(c) = num_groups * exponent * scl_base
result = num_groups * exponent * scl_base

# Print the final equation with each number included
print(f"The stable commutator length of c is calculated by the equation:")
print(f"{num_groups} * {exponent} * {scl_base} = {result}")
