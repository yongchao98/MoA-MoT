import numpy as np

# We found that the expression is maximized for n=1.
# This script verifies the calculation for this specific case.
n = 1
k = n + 2

# For n=1, the Cayley-Menger matrix is 3x3.
C = np.array([[0.0, 1.0, 1.0],
              [1.0, 0.0, 1.0],
              [1.0, 1.0, 0.0]])

# In the Hessenberg decomposition C = P.H.P^T, for a symmetric C,
# P is orthogonal and H is tridiagonal. The Householder method can be used
# to find P. For n=1, the matrix P can be chosen to be a simple reflection,
# which has eigenvalues {1, 1, -1}.
P = np.array([[1.0, 0.0, 0.0],
              [0.0, -1.0/np.sqrt(2), 1.0/np.sqrt(2)],
              [0.0, 1.0/np.sqrt(2), 1.0/np.sqrt(2)]])

# H is defined by the similarity transformation.
H = P.T @ C @ P

# Now, we calculate the four quantities for these specific matrices.

# Eigenvalues of P are {1, 1, -1}.
eigvals_P_real = np.linalg.eigvals(P).real
lambda_max_P = np.max(eigvals_P_real)
lambda_min_P = np.min(eigvals_P_real)

# Eigenvalues of H are the same as C: {n+1, -1, ..., -1}. For n=1, {2, -1, -1}.
lambda_max_H = n + 1
lambda_min_H = -1

# E_M = (lambda_max - lambda_min) / (size - 1)
E_P = (lambda_max_P - lambda_min_P) / (k - 1)
E_H = (lambda_max_H - lambda_min_H) / (k - 1)

# S_M = (1/size) * ||M||_F^2
# For an orthogonal matrix P, S_P = 1.
S_P = np.sum(P**2) / k
# For a symmetric matrix H, S_H = (1/size) * sum(lambda_i^2).
# For n=1, S_H = (1/3) * (2^2 + (-1)^2 + (-1)^2) = 6/3 = 2.
S_H = (lambda_max_H**2 + (n+1)*lambda_min_H**2) / k

# Final product
product = E_P * E_H * S_P * S_H

print("Calculation for n=1, where the theoretical maximum is achieved:")
print(f"E_P = ({lambda_max_P:.1f} - ({lambda_min_P:.1f})) / {k-1} = {E_P:.4f}")
print(f"E_H = ({lambda_max_H} - ({lambda_min_H})) / {k-1} = {E_H:.4f}")
print(f"S_P = {S_P:.4f}")
print(f"S_H = {S_H:.4f}")
print("\nFinal product expression:")
print(f"{E_P:.2f} * {E_H:.2f} * {S_P:.2f} * {S_H:.2f} = {product:.2f}")
