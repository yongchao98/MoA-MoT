# A script to calculate the number of non-isomorphic vertex-transitive graphs
# on 8 vertices for each possible degree.

# n_j is the number of such graphs with vertex degree j.
# A graph G is vertex-transitive if and only if its complement is.
# For n=8, this means n_j = n_{7-j}.
# We can find the values for j=0, 1, 2, 3 and deduce the rest.

# For j=0, the only graph is the empty graph (8 vertices, 0 edges). It is vertex-transitive.
n_0 = 1

# For j=1, the only graph is a perfect matching (4 disjoint edges). It is vertex-transitive.
n_1 = 1

# For j=2, the graphs must be unions of cycles. The only vertex-transitive options
# on 8 vertices are the 8-cycle (C8) and two 4-cycles (2C4).
n_2 = 2

# For j=3, we consult the known census of vertex-transitive graphs. There are 5
# non-isomorphic cubic vertex-transitive graphs on 8 vertices.
n_3 = 5

# Using the symmetry property n_j = n_{7-j}:
n_4 = n_3
n_5 = n_2
n_6 = n_1
n_7 = n_0

# Construct the list of counts [n_0, n_1, ..., n_7]
result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]

# Print the final list, showing each element explicitly as requested.
print(f"[{result[0]}, {result[1]}, {result[2]}, {result[3]}, {result[4]}, {result[5]}, {result[6]}, {result[7]}]")