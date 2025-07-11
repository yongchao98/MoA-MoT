# Define the parameters of the problem
num_groups = 19
exponent = 30

# The stable commutator length (scl) of a single basic commutator
# c_i = [a_i, b_i] in a free group F_i is 1/2.
scl_of_single_commutator = 0.5

# According to the homogeneity property of scl, scl(g^n) = n * scl(g).
# So, the scl of c_i^30 in the group F_i is:
scl_of_powered_commutator = exponent * scl_of_single_commutator

# According to the additivity property of scl over free products,
# the scl of the product c = c_1^30 * ... * c_19^30 is the sum
# of the individual scl values.
total_scl = num_groups * scl_of_powered_commutator

# Print the final equation with all the numbers involved in the calculation
print("The final calculation is based on the formula:")
print("Total scl = (Number of groups) * (Exponent) * (scl of a single commutator)")
print("\nFinal Equation:")
print(f"{num_groups} * {exponent} * {scl_of_single_commutator} = {total_scl}")