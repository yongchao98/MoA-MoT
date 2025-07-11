# Define the parameters of the barbell graph
total_nodes = 10
# A standard 10-node barbell graph consists of two cliques of 5 nodes each.
clique_size = 5

# Calculate the degree of the nodes connected by the bottleneck edge.
# Each of these two nodes is connected to (clique_size - 1) nodes in its own clique,
# plus 1 node in the other clique via the bridge.
deg_bottleneck_node = (clique_size - 1) + 1

# The bottleneck edge connects two such nodes. Let's call them u and v.
deg_u = deg_bottleneck_node
deg_v = deg_bottleneck_node

# The probability of sampling an edge (u, v) in a single step of uniform gossiping is:
# P = P(u is chosen) * P(u contacts v) + P(v is chosen) * P(v contacts u)
# P = (1 / total_nodes) * (1 / deg(u)) + (1 / total_nodes) * (1 / deg(v))
prob = (1 / total_nodes) * (1 / deg_u) + (1 / total_nodes) * (1 / deg_v)

# Print the result showing the full equation
print(f"A 10-node barbell graph has two K{clique_size} cliques.")
print(f"The degree of each node on the bottleneck edge is {deg_bottleneck_node}.")
print("The probability of sampling the bottleneck edge is calculated as:")
print(f"P = (1/{total_nodes}) * (1/{deg_u} + 1/{deg_v})")
print(f"P = {prob}")
