# Let 'd' represent the dimension of the qudit quantum system.
d_symbol = 'd'

# The maximal rank is d-squared.
# The only explicit number in this formula is the exponent, which is 2.
exponent = 2

# Print the final answer and the numbers in its equation.
print("For a Pauli channel acting on a qudit of dimension 'd',")
print("the maximal rank of its complementary channel is given by the equation:")
print(f"Maximal Rank = {d_symbol} ** {exponent}")
print("\n# The numbers in the final equation are:")
print(f"The base is the symbolic dimension: '{d_symbol}'")
print(f"The exponent is: {exponent}")