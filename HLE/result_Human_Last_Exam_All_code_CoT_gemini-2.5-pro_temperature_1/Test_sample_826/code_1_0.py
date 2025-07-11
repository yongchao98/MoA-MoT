import networkx as nx

# Step 1: Create two simple graphs, G = K2 and I = K1.
# K2 is a complete graph with 2 vertices (a single edge).
# K1 is a complete graph with 1 vertex (an isolated point).
# K1 is the only candidate for a multiplicative identity.
G = nx.complete_graph(2)
G.name = "G = K2"
I = nx.complete_graph(1)
I.name = "I = K1"

# Step 2: Compute the tensor product P = G ⊗ I.
P = nx.tensor_product(G, I)
P.name = "P = G ⊗ I"

# Step 3: Compare the original graph G with the product P.
# For I to be a multiplicative identity, P must be isomorphic to G.
# We can check this by comparing their number of vertices and edges.
print(f"Properties of the original graph {G.name}:")
print(f"Number of vertices: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")
print("-" * 30)
print(f"Properties of the product graph {P.name}:")
print(f"Number of vertices: {P.number_of_nodes()}")
print(f"Number of edges: {P.number_of_edges()}")
print("-" * 30)

# Step 4: Conclude whether they are isomorphic.
are_isomorphic = nx.is_isomorphic(G, P)

print(f"Is {G.name} isomorphic to {P.name}? {are_isomorphic}")
if not are_isomorphic:
    print("Since G is not isomorphic to G ⊗ I, I=K1 is not a multiplicative identity.")
    print("This shows that (G, ⊗) is not a monoid.")
