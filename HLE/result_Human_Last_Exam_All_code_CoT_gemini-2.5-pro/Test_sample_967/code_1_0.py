#
# Step 1: Define the properties of the Petersen graph
#
petersen_vertices = 10
petersen_edges = 15
petersen_degree = 3

#
# Step 2: Calculate the properties of the line graph of the Petersen graph (L(P))
# This is the underlying graph Gamma for our graph of groups.
#
# The number of vertices in L(P) is the number of edges in P.
num_vertices = petersen_edges

# The degree of each vertex in L(P) is (d-1) + (d-1) for a d-regular graph P.
degree_in_line_graph = (petersen_degree - 1) + (petersen_degree - 1)
# The number of edges in L(P) is (num_vertices * degree) / 2.
num_edges = (num_vertices * degree_in_line_graph) // 2

#
# Step 3: Define the l2-Betti numbers of the vertex and edge groups.
#
# As explained in the problem description analysis, the groups Ng are fundamental
# groups of finite covers of hyperbolic 3-manifolds. Their first l2-Betti number
# is 0.
# The vertex groups G_v are either Ng or free products of Ng's.
# The first l2-Betti number of a free product of infinite groups is the sum
# of their individual l2-Betti numbers.
# Therefore, the sum of the first l2-Betti numbers of all vertex groups is 0.
sum_betti1_vertex_groups = 0

# The edge groups are factors of the vertex groups, which are the Ng's themselves.
# Thus, their first l2-Betti numbers are also 0.
sum_betti1_edge_groups = 0

#
# Step 4: Apply the formula for the first l2-Betti number of a graph of groups.
# Formula: beta_1(G) = sum(beta_1(G_v)) - sum(beta_1(G_e)) + |E| - |V| + 1
#
betti_1_G = sum_betti1_vertex_groups - sum_betti1_edge_groups + num_edges - num_vertices + 1

#
# Step 5: Print the final equation and the result.
#
print("The first l2-Betti number is calculated using the formula:")
print("beta_1(G) = sum(beta_1(G_v)) - sum(beta_1(G_e)) + |E| - |V| + 1")
print("\nSubstituting the values:")
print(f"|V| (vertices in L(Petersen)) = {num_vertices}")
print(f"|E| (edges in L(Petersen)) = {num_edges}")
print(f"sum(beta_1(G_v)) = {sum_betti1_vertex_groups}")
print(f"sum(beta_1(G_e)) = {sum_betti1_edge_groups}")

print("\nThe final equation is:")
print(f"beta_1(G) = {sum_betti1_vertex_groups} - {sum_betti1_edge_groups} + {num_edges} - {num_vertices} + 1")
print(f"beta_1(G) = {betti_1_G}")
