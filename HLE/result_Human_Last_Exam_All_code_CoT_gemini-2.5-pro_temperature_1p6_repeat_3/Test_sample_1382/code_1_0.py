import numpy as np

# Define the graph K4
n = 4  # Number of nodes
# Edges of K4: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
m = len(edges)  # Number of edges

# --- Step 1: Construct the incidence matrix B ---
B = np.zeros((n, m))
for j, edge in enumerate(edges):
    u, v = edge
    B[u, j] = 1
    B[v, j] = 1

# --- Step 2: Calculate null(B^T B) ---
# null(B^T B) = m - rank(B^T B)
BtB = B.T @ B
rank_BtB = np.linalg.matrix_rank(BtB)
null_BtB = m - rank_BtB

print(f"Graph: K4")
print(f"Number of nodes (n): {n}")
print(f"Number of edges (m): {m}")
# For K4, it is connected (k=1) and not bipartite (z=0)
z = 0 
print(f"Number of bipartite components (z): {z}")

# We can also calculate it from the formula m - n + z
null_from_formula = m - n + z
print(f"Analytic nullity(B^T B) = m - n + z = {m} - {n} + {z} = {null_from_formula}")
print(f"Computed nullity(B^T B) from rank: {null_BtB}")
print("-" * 20)

# --- Step 3: Construct the Laplacian L and find its largest eigenvalue ---
A = np.zeros((n, n))
for edge in edges:
    u, v = edge
    A[u, v] = 1
    A[v, u] = 1

D = np.diag(np.sum(A, axis=1))
L = D - A

eigenvalues, _ = np.linalg.eigh(L)
lambda_n = eigenvalues[-1]

print("Laplacian Matrix L for K4:")
print(L)
print(f"\nEigenvalues of L: {np.round(eigenvalues, 4)}")
print(f"Largest eigenvalue λ_n(G): {lambda_n:.4f}")
print("-" * 20)

# --- Step 4: Verify the given relation ---
lhs = null_BtB
rhs = lambda_n / 2
print("Verifying the property: null(B^T B) = λ_n(G) / 2")
print(f"LHS = null(B^T B) = {lhs}")
print(f"RHS = λ_n(G) / 2 = {lambda_n:.4f} / 2 = {rhs:.4f}")
print(f"Is LHS == RHS? {np.isclose(lhs, rhs)}")
print("-" * 20)

# --- Step 5: Check the answer choices for K4 ---
print("Checking Answer Choices for K4:")
# C. The graph is connected.
# K4 is connected, so k=1.
k = 1
print(f"Choice C: The graph is connected (k=1). This is TRUE for K4.")
# D. The graph has exactly lambda_n(G)/2 connected components.
# k == lambda_n/2  => 1 == 4/2 = 2.
print(f"Choice D: k = λ_n/2?  => {k} == {rhs}? {k == rhs}. This is FALSE for K4.")
