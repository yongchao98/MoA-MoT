import sympy

# Let d be the degree of the graph (number of neighbors for any vertex)
# For a 2D torus, d = 4
d = sympy.Integer(4)

# Let k be the number of common neighbors between vertex 0 and vertex x_0
# As determined in the analysis, x_0 is a diagonal neighbor, so k = 2
k = sympy.Integer(2)

# The limiting conditional probability is given by the ratio of "private" neighbors
# to the total number of neighbors, which is (d-k)/d.
result_fraction = (d - k) / d
result_float = float(result_fraction)

print(f"The degree of each vertex is d = {d}.")
print(f"The number of common neighbors is k = {k}.")
print("The limit of the conditional probability is (d - k) / d.")
print(f"P = ({d} - {k}) / {d} = {result_fraction}")
print(f"The numerical value is {result_float}.")
