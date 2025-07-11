# The number of vertices in the line graph of the Petersen graph.
num_vertices = 15
# The number of edges in the line graph of the Petersen graph.
num_edges = 30

# The first Betti number of the underlying graph Y, which contributes to the final result.
# beta_1(Y) = |E| - |V| + 1
beta_1_graph = num_edges - num_vertices + 1

# Calculate the sum of the first L2-Betti numbers of the vertex groups.
# For v_1, the group is N_100, and its Betti number is 0.
sum_beta_1_vertex_groups = 0
# For v_i where i is from 2 to 15, the group is a free product of i-1 infinite groups
# with vanishing first L2-Betti numbers.
# The formula for the first L2-Betti number of a free product H_1 * ... * H_k is:
# sum(beta_1(H_j)) + k - 1
# In our case, this simplifies to (i-1) - 1 = i - 2.
for i in range(2, num_vertices + 1):
    sum_beta_1_vertex_groups += (i - 2)

# Calculate the sum of the first L2-Betti numbers of the edge groups.
# Each edge group is some N_g, for which the first L2-Betti number is 0.
# So the total sum over all edges is 0.
sum_beta_1_edge_groups = 0

# The formula for the first l2-betti number of the fundamental group G of a graph of groups is:
# beta_1(G) = sum(beta_1(G_v)) - sum(beta_1(G_e)) + beta_1(Y)
final_betti_number = sum_beta_1_vertex_groups - sum_beta_1_edge_groups + beta_1_graph

# Print the final equation with all its components.
print("The first l2-Betti number is computed from the formula:")
print("Sum(beta_1(G_v)) - Sum(beta_1(G_e)) + beta_1(Y)")
print(f"The equation is: {sum_beta_1_vertex_groups} - {sum_beta_1_edge_groups} + {beta_1_graph} = {final_betti_number}")
