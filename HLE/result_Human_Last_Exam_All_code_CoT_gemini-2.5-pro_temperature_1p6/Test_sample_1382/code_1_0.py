import networkx as nx
import numpy as np

# Step 1: Create the graph K4 (complete graph on 4 nodes)
G = nx.complete_graph(4)

# Step 2: Get the graph's properties
n = G.number_of_nodes()
m = G.number_of_edges()
k = nx.number_connected_components(G)

# Determine kb, the number of bipartite connected components
k_b = 0
for component_nodes in nx.connected_components(G):
    subgraph = G.subgraph(component_nodes)
    if nx.is_bipartite(subgraph):
        k_b += 1

# Step 3: Calculate the terms of the given equation

# LHS: nullity(B^T B) = m - n + kb
nullity_b_transpose_b = m - n + k_b

# RHS: lambda_n(L) / 2
# Get the graph Laplacian L
L = nx.laplacian_matrix(G).toarray()
# Calculate its eigenvalues. Use eigvalsh since L is symmetric.
eigenvalues = np.linalg.eigvalsh(L)
# Find the largest eigenvalue lambda_n
lambda_n = np.max(eigenvalues)
half_lambda_n = lambda_n / 2

# Step 4: Verify that K4 satisfies the condition and output the equation
print("Analyzing the graph K4:")
print(f"Number of nodes (n): {n}")
print(f"Number of edges (m): {m}")
print(f"Number of bipartite components (k_b): {k_b}")
print(f"Largest Laplacian eigenvalue (lambda_n): {lambda_n:.4f}")
print("-" * 30)

print("Checking if the condition null(B^T B) = lambda_n(G)/2 holds.")
print(f"The equation is: {m} - {n} + {k_b} = {lambda_n:.4f} / 2")
print(f"LHS = {nullity_b_transpose_b}")
print(f"RHS = {half_lambda_n}")
print(f"Does the condition hold? {np.isclose(nullity_b_transpose_b, half_lambda_n)}")
print("-" * 30)


# Step 5: Test the answer choices for K4. Let's test option D explicitly.
# D. The graph has exactly lambda_n(G)/2 connected components
print("Testing Answer Choice D:")
num_components = k
print(f"Number of connected components (k): {num_components}")
print(f"Value of lambda_n(G)/2: {half_lambda_n}")
print(f"Is k == lambda_n(G)/2? {num_components == half_lambda_n}")

# Analysis of other options for K4:
# A. If you drop lambda_n/2 = 2 edges, there will be at least two nodes with degree <= 1.
#    If we remove two non-adjacent edges from K4, we get a 4-cycle (C4),
#    where every node has degree 2. So, choice A is false.
# B. The graph has at least lambda_n/2 connected components (k >= lambda_n/2).
#    For K4, this is 1 >= 2, which is false.
# C. The graph is connected. K4 is connected, so this is true for K4.
#    However, a graph of 4 isolated nodes also satisfies the initial condition
#    (LHS=0-4+4=0, RHS=0/2=0), but is not connected. So C is not always true.
