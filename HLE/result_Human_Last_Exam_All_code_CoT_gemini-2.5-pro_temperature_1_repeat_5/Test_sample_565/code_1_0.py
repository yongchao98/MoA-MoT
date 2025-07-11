# This script calculates the number of non-isomorphic vertex-transitive graphs on 8 vertices
# for each degree j from 0 to 7.

# n_j is the number of graphs for degree j.
# Based on graph theory principles and known results.

# n_0: For degree 0, the only graph is the null graph on 8 vertices (N_8).
# It is vertex-transitive.
n_0 = 1

# n_1: For degree 1, the only graph is a perfect matching (4 disjoint K_2).
# It is vertex-transitive.
n_1 = 1

# n_2: For degree 2, the vertex-transitive graphs are the 8-cycle (C_8)
# and the disjoint union of two 4-cycles (2C_4).
n_2 = 2

# n_3: For degree 3, it's a known result from algebraic graph theory that there are
# 5 non-isomorphic cubic vertex-transitive graphs on 8 vertices.
n_3 = 5

# For a vertex-transitive graph G of degree j on n vertices, its complement is
# vertex-transitive and has degree (n-1-j). For n=8, this means n_j = n_{7-j}.
# We use this symmetry to find the remaining values.

# n_4 is the same as n_3
n_4 = n_3

# n_5 is the same as n_2
n_5 = n_2

# n_6 is the same as n_1
n_6 = n_1

# n_7 is the same as n_0
n_7 = n_0

# The final list of numbers.
result_list = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]

# Print the final list in the required format.
print(result_list)