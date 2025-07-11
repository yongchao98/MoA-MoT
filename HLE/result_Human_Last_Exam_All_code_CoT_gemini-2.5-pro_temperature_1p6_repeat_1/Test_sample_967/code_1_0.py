# Step 1: Compute properties of the underlying graph Gamma = L(Petersen)
# The Petersen graph has 10 vertices, each with degree 3.
num_vertices_petersen = 10
degree_petersen = 3

# The number of vertices in the line graph is the number of edges in the original graph.
# Petersen graph has (10 * 3) / 2 = 15 edges.
num_vertices_gamma = 15

# The number of edges in the line graph is sum over vertices v of (d(v) choose 2).
# For the 3-regular Petersen graph, this is 10 * (3 choose 2) = 10 * 3 = 30.
num_edges_gamma = 30

# The first Betti number of the graph Gamma (since it's connected).
# This is also the first l2-Betti number of Gamma.
b1_gamma = num_edges_gamma - num_vertices_gamma + 1

# Step 2: Compute the sum of the first l2-Betti numbers of the vertex groups.
# For g >= 2, beta_1^{(2)}(N_g) = 0.
# The vertex group G_v1 is N_100, so its b1 is 0.
b1_G_v1 = 0

# The vertex groups G_vi for i=2..15 are free products of N_g's.
# beta_1^{(2)}(G_vi) = i - 2.
# We sum (i-2) for i from 2 to 15. This is equivalent to summing j from 0 to 13.
sum_b1_vertex_groups = b1_G_v1 + sum(i - 2 for i in range(2, 16))

# Step 3: Compute the sum of the first l2-Betti numbers of the edge groups.
# Each edge group G_e is some N_g, so its b1 is 0.
# There are 30 edges.
sum_b1_edge_groups = 30 * 0

# Step 4: Combine the results using the main formula.
final_betti_number = sum_b1_vertex_groups - sum_b1_edge_groups + b1_gamma

# Print the final equation with each number.
print(f"{sum_b1_vertex_groups} - {sum_b1_edge_groups} + {b1_gamma} = {final_betti_number}")
print("The first l2-betti number is:")
print(final_betti_number)