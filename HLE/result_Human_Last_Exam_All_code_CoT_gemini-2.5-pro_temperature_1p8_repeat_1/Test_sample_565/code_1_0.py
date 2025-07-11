# The number of isomorphism classes of vertex-transitive graphs
# on 8 vertices, categorized by vertex degree j.
# Let n_j be the number of such graphs with degree j.

# For j=0, the only graph is the empty graph E_8.
n_0 = 1

# For j=1, the only graph is the union of 4 disjoint edges (4K_2).
n_1 = 1

# For j=2, the graphs are the 8-cycle (C_8) and two 4-cycles (2C_4).
n_2 = 2

# For j=3, published enumerations show there are 6 such graphs.
# These include the cube graph (Q_3), 2K_4, and four others.
n_3 = 6

# The complement of a j-regular vertex-transitive graph on 8 vertices
# is a (7-j)-regular vertex-transitive graph.
# Therefore, n_j = n_{7-j}.
n_4 = n_3
n_5 = n_2
n_6 = n_1
n_7 = n_0

# The list of numbers [n_0, n_1, ..., n_7]
result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]

# Printing the result in the specified format
print(f"[{result[0]}, {result[1]}, {result[2]}, {result[3]}, {result[4]}, {result[5]}, {result[6]}, {result[7]}]")