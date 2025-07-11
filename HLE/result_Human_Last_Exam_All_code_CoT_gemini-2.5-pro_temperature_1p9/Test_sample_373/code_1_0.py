# Define the parameters of the barbell graph and the gossiping process
total_nodes = 10
nodes_in_clique = 5

# Calculate the degree of the two nodes connected by the bottleneck edge.
# Each node is connected to (nodes_in_clique - 1) nodes in its own clique,
# plus 1 connection for the bottleneck edge itself.
degree_of_bottleneck_node_u = (nodes_in_clique - 1) + 1
degree_of_bottleneck_node_v = (nodes_in_clique - 1) + 1

# The probability P(u,v) of sampling the edge (u, v) is given by:
# P(u,v) = (1/N) * (1/d(u) + 1/d(v))
# where N is the total number of nodes, and d(u) and d(v) are the degrees.
probability = (1 / total_nodes) * (1 / degree_of_bottleneck_node_u + 1 / degree_of_bottleneck_node_v)

# Print the final equation with the numbers plugged in
print("The probability of sampling the bottleneck edge (u,v) is calculated with the formula:")
print("P(u,v) = (1 / Total Nodes) * (1 / degree(u) + 1 / degree(v))")
print("\nHere is the calculation with the specific numbers for this graph:")
print(f"(1 / {total_nodes}) * (1 / {degree_of_bottleneck_node_u} + 1 / {degree_of_bottleneck_node_v}) = {probability}")