# Define the parameters from the problem statement.
num_groups = 19
exponent = 30

# This is a known mathematical result for the stable commutator length (scl)
# of a single commutator [a, b] in the free group F_2 on generators {a, b}.
scl_of_single_commutator = 0.5

# The final calculation is based on the properties of scl:
# scl(c) = scl(product of c_i^30)
#        = sum(scl(c_i^30))          (by additivity over free products)
#        = sum(exponent * scl(c_i))  (by homogeneity)
#        = num_groups * exponent * scl(c_i)
total_scl = num_groups * exponent * scl_of_single_commutator

# Print the final equation with each number and the result.
# The result should be an integer, so we cast it.
print(f"The stable commutator length of c is given by the equation:")
print(f"{num_groups} * {exponent} * {scl_of_single_commutator} = {int(total_scl)}")