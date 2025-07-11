import math

# Step 1: Define the graph parameters
total_nodes = 10
num_cliques = 2

# For a standard, symmetric barbell graph, the nodes are split equally between the cliques.
nodes_per_clique = total_nodes // num_cliques

# Step 2: Calculate the number of edges in each component
# The number of edges in a complete graph (clique) with 'm' nodes is m * (m - 1) / 2
edges_in_clique1 = nodes_per_clique * (nodes_per_clique - 1) // 2
edges_in_clique2 = edges_in_clique1  # The cliques are identical

# The two cliques are connected by a single bottleneck edge.
bottleneck_edges = 1

# Step 3: Calculate the total number of edges in the graph
total_edges = edges_in_clique1 + edges_in_clique2 + bottleneck_edges

# Step 4: Calculate the probability
# In uniform sampling, the probability of picking any single edge is 1 / (total number of edges).
numerator = bottleneck_edges
denominator = total_edges
probability = numerator / denominator

# Final Output: Print the explanation and the final equation.
print(f"A barbell graph with {total_nodes} nodes consists of two {nodes_per_clique}-node complete graphs connected by a bridge.")
print(f"Number of edges in the first clique: {edges_in_clique1}")
print(f"Number of edges in the second clique: {edges_in_clique2}")
print(f"Number of bottleneck edges: {bottleneck_edges}\n")
print("The probability of sampling the bottleneck edge is the number of bottleneck edges divided by the total number of edges.")
print("The final equation is:")
print(f"P(bottleneck) = {numerator} / ({edges_in_clique1} + {edges_in_clique2} + {bottleneck_edges}) = {numerator} / {denominator}")
print(f"\nThe calculated probability is: {probability}")