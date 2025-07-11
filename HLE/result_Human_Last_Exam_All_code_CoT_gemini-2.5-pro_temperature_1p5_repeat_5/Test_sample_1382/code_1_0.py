import numpy as np

# Let's verify the reasoning with the 'two triangles sharing a vertex' graph.
# n = 5 nodes (0, 1, 2, 3, 4)
# m = 6 edges
# Triangle 1 on (0, 1, 2), Triangle 2 on (0, 3, 4). Node 0 is the shared vertex.
# Adjacency Matrix A
A = np.array([
    [0, 1, 1, 1, 1],
    [1, 0, 1, 0, 0],
    [1, 1, 0, 0, 0],
    [1, 0, 0, 0, 1],
    [1, 0, 0, 1, 0]
])

# Degree Matrix D
D = np.diag(np.sum(A, axis=1))

# Laplacian Matrix L = D - A
L = D - A

# Eigenvalues of L
eigenvalues = np.linalg.eigvalsh(L)
lambda_n = np.max(eigenvalues)

# Graph properties
n = 5
m = 6
c = 1

# Check the condition from the problem statement
# LHS: cyclomatic number
cyclomatic_number = m - n + c
# RHS: lambda_n / 2
lambda_n_div_2 = lambda_n / 2

print("Graph: Two triangles sharing a vertex")
print(f"n = {n}, m = {m}, c = {c}")
print(f"Maximum degree d_max = {int(np.max(np.diag(D)))}")
print(f"LHS of the equation (m - n + c) = {cyclomatic_number}")
print(f"Largest eigenvalue lambda_n = {lambda_n:.4f}")
print(f"RHS of the equation (lambda_n / 2) = {lambda_n_div_2:.4f}")
print("\nConclusion from the numeric example:")
print(f"The equation requires {cyclomatic_number} = {lambda_n_div_2:.4f}, which is {cyclomatic_number == lambda_n_div_2}.")
print("More importantly, we test the bound lambda_n >= d_max + 1.")
print(f"Is {lambda_n:.4f} >= {int(np.max(np.diag(D)))} + 1? This is {lambda_n >= np.max(np.diag(D)) + 1}.")
print("The property requires lambda_n = 2 * (m-n+c) = 2 * 2 = 4.")
print("But the spectral bound requires lambda_n >= 5.")
print("This is a contradiction, so this graph (with m>n) cannot have the property.")
print("This supports the conclusion that having the property implies m <= n, which is statement B.")
