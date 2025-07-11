import numpy as np

# Let's define a matrix A, and vectors b, x
A = np.array([[1.0, 0.0],
              [0.0, 2.0],
              [3.0, 1.0]])

x = np.array([[2.0],
              [1.0]])

b = np.array([[5.0],
              [5.0],
              [5.0]])

# We are looking for the minimum norm matrix E such that x solves the
# least-squares problem for (A+E)x = b.
# A simplified, but valid, special case is when (A+E)x = b exactly.
# This leads to Ex = b - Ax. Let r = b - Ax.

# Calculate Ax
Ax = A @ x

# Calculate the residual r = b - Ax
r = b - Ax

# The minimum Frobenius norm E satisfying Ex = r is given by
# E = r * x^T / (x^T * x)
x_T = x.T
xTx = x_T @ x
scalar_denominator = xTx[0, 0]

# Calculate E
E = (r @ x_T) / scalar_denominator

# Calculate the rank of E
rank_E = np.linalg.matrix_rank(E)

print("Given:")
print("A =\n", A)
print("\nx =\n", x)
print("\nb =\n", b)
print("\nWe consider the case where (A+E)x = b, which gives Ex = b - Ax.")
print("First, we calculate r = b - Ax:")
print("r =\n", r)
print("\nThe minimal Frobenius norm E is given by the equation: E = r * x^T / (x^T*x)")
print("Let's output each number in the final equation:")
print("E =\n", E)
print("\nr =\n", r)
print("\nx^T =\n", x_T)
print("\nx^T*x = ", scalar_denominator)

# Verification
print("\nDoes this E satisfy Ex = r?")
print("Ex =\n", E @ x)
print("r =\n", r)


print(f"\nThe rank of the resulting matrix E is: {rank_E}")
