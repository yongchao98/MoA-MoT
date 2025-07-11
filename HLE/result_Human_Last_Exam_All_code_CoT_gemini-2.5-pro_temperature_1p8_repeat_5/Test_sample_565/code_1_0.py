# The task is to calculate the numbers n_j of isomorphism classes of
# vertex-transitive graphs with 8 vertices and vertex degree j for j = 0, ..., 7.

# Initialize a list to store the counts n_j for j = 0 to 7.
n_counts = [0] * 8

# For a graph on n vertices, vertex-transitivity implies the graph is regular.
# For n=8, we consider degrees j from 0 to 7.

# n_0: Number of 0-regular vertex-transitive graphs on 8 vertices.
# The only such graph is the empty graph (8 isolated vertices), which is vertex-transitive.
n_counts[0] = 1

# n_1: Number of 1-regular vertex-transitive graphs on 8 vertices.
# Any 1-regular graph on 8 vertices is a perfect matching (4K_2).
# It is unique up to isomorphism and is vertex-transitive.
n_counts[1] = 1

# n_2: Number of 2-regular vertex-transitive graphs on 8 vertices.
# 2-regular graphs are disjoint unions of cycles. For 8 vertices, the
# vertex-transitive options are a single 8-cycle (C_8) or two 4-cycles (2C_4).
# This gives two non-isomorphic graphs.
n_counts[2] = 2

# n_3: Number of 3-regular (cubic) vertex-transitive graphs on 8 vertices.
# This is a known result from graph enumeration. There are 5 such non-isomorphic graphs.
# These include one disconnected graph (2K_4) and four connected graphs (one being the cube, Q_3).
n_counts[3] = 5

# The number of j-regular vertex-transitive graphs on n vertices is equal to the number of
# (n-1-j)-regular vertex-transitive graphs, due to graph complementation.
# For n=8, this means n_j = n_{7-j}. We use this to find the remaining values.

# n_4 = n_{7-4} = n_3
n_counts[4] = n_counts[3]

# n_5 = n_{7-5} = n_2
n_counts[5] = n_counts[2]

# n_6 = n_{7-6} = n_1
n_counts[6] = n_counts[1]

# n_7 = n_{7-7} = n_0
n_counts[7] = n_counts[0]

# Print the final list in the required format.
print(n_counts)