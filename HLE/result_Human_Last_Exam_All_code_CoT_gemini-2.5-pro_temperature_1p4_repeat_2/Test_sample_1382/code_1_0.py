import networkx as nx
import numpy as np

# Create a graph that satisfies the condition.
# A disjoint union of two 4-cycles is a good example.
G1 = nx.cycle_graph(4)
G2 = nx.cycle_graph(4)
G = nx.disjoint_union(G1, G2)

# 1. Calculate the properties for the cyclomatic number
n = G.number_of_nodes()
m = G.number_of_edges()
c = nx.number_connected_components(G)

# The cyclomatic number mu is null(B^T B)
mu = m - n + c

# 2. Calculate the largest eigenvalue of the Laplacian
L = nx.laplacian_matrix(G).toarray()
# Use eigvalsh for symmetric matrices like the Laplacian
eigenvalues = np.linalg.eigvalsh(L)
lambda_n = np.max(eigenvalues)

# The problem states that mu = lambda_n / 2. Let's verify.
print("Verifying the given condition: null(B^T B) = lambda_n(G) / 2")
print(f"Number of nodes (n): {n}")
print(f"Number of edges (m): {m}")
print(f"Number of connected components (c): {c}")
print(f"Cyclomatic number (mu = m - n + c): {mu}")
print(f"Largest Laplacian eigenvalue (lambda_n): {lambda_n:.4f}")
print(f"Equation from problem statement: {mu} = {lambda_n:.0f} / 2")
print("-" * 30)

# 3. What does this imply? We concluded it implies Choice D.
# Let's check Choice D: c = lambda_n / 2
print("Checking the implication from Choice D: c = lambda_n / 2")
# We use integers for the final check as the values are exact.
c_val = c
lambda_n_val = round(lambda_n)
val_from_D = lambda_n_val / 2

print(f"The number of components is {c_val}.")
print(f"The value of lambda_n / 2 is {lambda_n_val} / 2 = {val_from_D}.")

if c_val == val_from_D:
    print("\nConclusion: The statement from Choice D holds true for this graph.")
    print(f"The final equation is: {c_val} = {lambda_n_val} / 2")
else:
    print("\nConclusion: The statement from Choice D does NOT hold true.")
