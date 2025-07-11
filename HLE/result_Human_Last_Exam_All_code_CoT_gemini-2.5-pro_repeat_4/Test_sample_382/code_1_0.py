import numpy as np

# Let's define some example matrix A, vector b, and non-zero vector x.
# We choose them such that x is not already a solution for Ax = b.
A = np.array([
    [1.0, 2.0],
    [3.0, 4.0],
    [5.0, 6.0]
])

b = np.array([
    [7.0],
    [8.0],
    [9.0]
])

x = np.array([
    [1.0],
    [1.0]
])

# Step 1: Calculate the residual r = b - Ax
# This is the vector that Ex must equal.
Ax = A @ x
r = b - Ax

print("The equation for the minimal perturbation E is Ex = r")
print("where r = b - Ax\n")

# Step 2: The solution for E is r * x_pseudoinverse, which is (r @ x.T) / (x.T @ x)
# We will compute the components of this equation.
xtx = (x.T @ x).item() # The scalar value of x.T * x
r_xT = r @ x.T         # The outer product r * x.T

print(f"Vector r:\n{r}\n")
print(f"Vector x transpose:\n{x.T}\n")
print(f"Scalar x.T @ x: {xtx:.1f}\n")

# Calculate the matrix E with the minimum Frobenius norm
E = r_xT / xtx

print(f"Resulting matrix E:\n{E}\n")

# Step 3: Calculate the rank of E. Since r is non-zero and x is non-zero,
# the rank of their outer product should be 1.
rank_E = np.linalg.matrix_rank(E)

print(f"The rank of matrix E is: {rank_E}")

# We can also verify that (A+E)x = b
perturbed_Ax = (A + E) @ x
print(f"\nVerification: (A+E)x =\n{perturbed_Ax}")
print(f"Is (A+E)x equal to b? {np.allclose(perturbed_Ax, b)}")