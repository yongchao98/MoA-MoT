import numpy as np

# Define the graph K4
n = 4  # Number of nodes
m = 6  # Number of edges

# Unoriented incidence matrix B for K4
# Edges are (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
B = np.array([
    [1, 1, 1, 0, 0, 0],  # Node 0
    [1, 0, 0, 1, 1, 0],  # Node 1
    [0, 1, 0, 1, 0, 1],  # Node 2
    [0, 0, 1, 0, 1, 1]   # Node 3
])

# Adjacency matrix A for K4
A = np.array([
    [0, 1, 1, 1],
    [1, 0, 1, 1],
    [1, 1, 0, 1],
    [1, 1, 1, 0]
])

# Degree matrix D for K4
D = np.diag(np.sum(A, axis=1))

# Laplacian matrix L = D - A
L = D - A

# --- Calculate the Left-Hand Side (LHS) ---
# null(B^T B) = m - n + c_b
# For K4, c=1 and it's not bipartite, so c_b=0
c_b = 0
lhs = m - n + c_b

# --- Calculate the Right-Hand Side (RHS) ---
# Calculate eigenvalues of L and find the largest one
eigenvalues, _ = np.linalg.eigh(L)
lambda_n = np.max(eigenvalues)
rhs = lambda_n / 2

# --- Final Output ---
print(f"Graph: K4 (Complete graph on 4 vertices)")
print(f"n (nodes) = {n}")
print(f"m (edges) = {m}")
print(f"c_b (bipartite connected components) = {c_b}")
print("-" * 20)
print("Equation: null(B^T B) = lambda_n(G) / 2")
print("Equivalent to: m - n + c_b = lambda_n(G) / 2")
print("-" * 20)
print(f"LHS calculation: m - n + c_b = {m} - {n} + {c_b} = {lhs}")
print(f"RHS calculation: lambda_n(G) / 2 = {lambda_n:.2f} / 2 = {rhs:.2f}")

if np.isclose(lhs, rhs):
    print("\nThe equation holds true for K4.")
else:
    print("\nThe equation does not hold for K4.")

print("\nNow we check the answer choices for K4:")
c = 1 # K4 is connected
print(f"B: Is c >= lambda_n/2? Is {c} >= {rhs}? {'True' if c >= rhs else 'False'}")
print(f"C: Is the graph connected? {'True' if c == 1 else 'False'}")
print(f"D: Is c = lambda_n/2? Is {c} == {rhs}? {'True' if c == rhs else 'False'}")
print("\nSince options B and D are False for this example, they cannot be the correct answer.")
print("A separate example (K4 U K1) shows that C is also not always true.")
print("Therefore, none of the options A, B, C, or D is a necessary consequence.")
