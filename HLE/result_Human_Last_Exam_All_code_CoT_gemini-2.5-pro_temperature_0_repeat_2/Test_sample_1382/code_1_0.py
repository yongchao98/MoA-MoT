import networkx as nx
import numpy as np

# 1. Create a graph that satisfies the property.
# A graph made of two disjoint 4-cycles (C4) will work.
G = nx.Graph()
# Add the first C4 on nodes 0, 1, 2, 3
G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0)])
# Add the second C4 on nodes 4, 5, 6, 7
G.add_edges_from([(4, 5), (5, 6), (6, 7), (7, 4)])

# 2. Calculate the number of connected components, c.
c = nx.number_connected_components(G)

# 3. Calculate the graph Laplacian L.
L = nx.laplacian_matrix(G).toarray()

# 4. Calculate the eigenvalues of L.
# We use eigvalsh because the Laplacian is a symmetric matrix.
eigenvalues = np.linalg.eigvalsh(L)

# 5. Find the largest eigenvalue, lambda_n.
lambda_n = eigenvalues[-1]

# 6. Print the results and check the relationship.
print(f"The graph has n = {G.number_of_nodes()} nodes and m = {G.number_of_edges()} edges.")
print(f"Number of connected components (c): {c}")
print(f"Largest Laplacian eigenvalue (lambda_n): {lambda_n:.4f}")
print(f"Value of lambda_n / 2: {lambda_n / 2:.4f}")
print("\nChecking the condition from Answer D: c = lambda_n / 2")
print(f"Equation: {c} = {lambda_n} / 2")
print(f"Result: {c} = {lambda_n / 2}")
is_satisfied = np.isclose(c, lambda_n / 2)
print(f"Is the condition satisfied? {is_satisfied}")
