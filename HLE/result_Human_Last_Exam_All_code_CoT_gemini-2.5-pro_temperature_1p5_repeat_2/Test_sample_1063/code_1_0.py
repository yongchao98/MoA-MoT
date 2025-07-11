import networkx as nx

# Define the actors (nodes)
actors = [
    "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
    "Charlie Cox", "Devery Jacobs", "Thomas Elms"
]

# Define the co-starring relationships (edges) based on the criteria
# 1. Aaron Ashmore & Emilia Jones in 'Locke & Key' (S1: 2020)
# 2. Krysten Ritter & Charlie Cox in 'The Defenders' (Miniseries: 2017)
# 3. Devery Jacobs & Thomas Elms in 'The Order' (S1: 2019)
edges = [
    ("Aaron Ashmore", "Emilia Jones"),
    ("Krysten Ritter", "Charlie Cox"),
    ("Devery Jacobs", "Thomas Elms")
]

# Create the graph
G = nx.Graph()
G.add_nodes_from(actors)
G.add_edges_from(edges)

# Analyze the graph properties
is_connected = nx.is_connected(G)
# A graph is acyclic if it is a forest.
is_acyclic = nx.is_forest(G)

# Print the findings
print(f"Nodes: {list(G.nodes())}")
print(f"Edges: {list(G.edges())}")
print(f"Number of connected components: {nx.number_connected_components(G)}")
print("-" * 20)
print(f"Is the graph connected? {is_connected}")
print(f"Is the graph acyclic? {is_acyclic}")
print("-" * 20)

# Determine the final answer based on the analysis
if not is_connected and is_acyclic:
    result_description = "A. Disconnected and acyclic"
elif not is_connected and not is_acyclic:
    result_description = "B. Disconnected and cyclic"
elif is_connected and is_acyclic:
    result_description = "C. Connected and acyclic"
elif is_connected and not is_acyclic:
    if nx.is_cycle(G, G.nodes()):
         result_description = "E. A cycle graph"
    else:
         result_description = "D. Connected and cyclic, but not a cycle graph"
else:
    result_description = "Could not determine graph type."

print(f"The graph is correctly described as: {result_description}")
print("<<<A>>>")