# This script calculates the number of isomorphism classes of vertex-transitive
# graphs on 8 vertices, categorized by vertex degree j.
# Let n_j be the number of such graphs with degree j.

# Step 1: Establish the symmetry property.
# A graph G on n=8 vertices is vertex-transitive and k-regular if and only if
# its complement is vertex-transitive and (7-k)-regular.
# This implies that n_k must be equal to n_{7-k}.

# Step 2: Manually enumerate the graphs for degrees k = 0, 1, 2, 3.

# For k=0: The only graph is the empty graph N_8. It is vertex-transitive.
n_0 = 1

# For k=1: The only graph is the union of 4 disjoint edges (4K_2). It is vertex-transitive.
n_1 = 1

# For k=2: The graphs must be unions of cycles. The vertex-transitive options are:
# 1. The 8-cycle (C_8).
# 2. Two disjoint 4-cycles (2C_4).
# These are non-isomorphic.
n_2 = 2

# For k=3: The vertex-transitive options are:
# 1. The connected cubical graph (Q_3).
# 2. The disconnected graph made of two K_4 components (2K_4).
# These are non-isomorphic.
n_3 = 2

# Step 3: Apply the symmetry property n_k = n_{7-k} to find the remaining counts.
n_4 = n_3
n_5 = n_2
n_6 = n_1
n_7 = n_0

# Step 4: Assemble the final list.
# The list contains the number of graphs for each degree j from 0 to 7.
result_list = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]

# Print the final list as the answer.
print(result_list)