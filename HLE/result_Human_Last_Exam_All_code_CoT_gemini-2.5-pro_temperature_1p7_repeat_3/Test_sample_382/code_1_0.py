import numpy as np

# Set up the problem with a matrix A, and vectors b and x.
# We choose them such that x is not a solution to Ax = b.
A = np.array([[1, 0],
              [2, 1],
              [0, 1]])

b = np.array([[5],
              [5],
              [5]])

x = np.array([[1],
              [2]])

print("--- Given Matrices and Vectors ---")
print("Matrix A:\n", A)
print("\nVector b:\n", b)
print("\nNon-zero vector x:\n", x)

# --- Step 1: Calculate the residual r = b - Ax ---
# This is the vector that Ex must equal.
Ax = A @ x
r = b - Ax

print("\n--- Step-by-Step Calculation of E ---")
print("The goal is to find E with minimum ||E||_F such that (A+E)x = b")
print("This simplifies to finding E such that Ex = b - Ax = r")
print("\n1. Calculate the residual r = b - Ax:")
print("A @ x =\n", Ax)
print("\nr = b - (A @ x) =\n", r)

# --- Step 2: Use the formula E = r * x^T / (x^T * x) ---
# The solution E is the outer product of r and x, scaled by the squared norm of x.
xt = x.T
xtx = xt @ x
rxt = r @ xt

E = rxt / xtx

print("\n2. The solution is given by the formula E = (r * x^T) / (x^T * x)")
print("\nx^T =\n", xt)
print("\nx^T * x =\n", xtx)
print("\nr * x^T =\n", rxt)
print("\nE = (r * x^T) / (x^T * x) =\n", E)

# --- Step 3: Calculate the rank of E ---
rank_E = np.linalg.matrix_rank(E)

print("\n--- Resulting Rank ---")
print(f"The rank of the resulting matrix E is: {rank_E}")

# --- Step 4: Verification ---
# Let's check if (A+E)x = b is satisfied.
A_plus_E = A + E
result = A_plus_E @ x

print("\n--- Verification ---")
print("Let's verify that with our calculated E, the equation (A+E)x = b holds true.")
print("A + E =\n", A_plus_E)
print("\n(A + E) @ x =\n", result)
print("\nThis matches our target vector b.")

print("\nAs shown by the derivation and the example, the rank of E is 1 when r is non-zero.")
print("If r were zero (meaning Ax=b initially), E would be the zero matrix with rank 0.")
print("Therefore, the greatest possible rank of E is 1.")
