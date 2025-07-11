import numpy as np
from scipy.linalg import expm

def skew(v):
    """
    Creates a skew-symmetric matrix from a 3-element vector.
    v -> [[0, -v_z, v_y], [v_z, 0, -v_x], [-v_y, v_x, 0]]
    """
    if len(v) != 3:
        raise ValueError("Input vector must have 3 elements.")
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

# Step 1: Define initial state and covariance before the reset.
# The attitude deviation vector estimate, delta.
delta_vec = np.array([0.1, -0.05, 0.02])
# The corresponding 3x3 covariance matrix, Sigma.
Sigma = np.array([[0.01, 0.001, -0.002],
                  [0.001, 0.009, 0.003],
                  [-0.002, 0.003, 0.012]])

# Step 2: Compute the matrices required for the transformation equation.
# The equation is: Sigma' = exp(-Delta) * Sigma * exp(Delta)

# Compute the skew-symmetric matrix, Delta, from the deviation vector.
Delta = skew(delta_vec)

# Compute the transformation matrix G = exp(-Delta) using the matrix exponential.
G = expm(-Delta)

# Compute the transpose of the transformation matrix, G_T = exp(Delta).
# This is valid because (exp(A))^T = exp(A^T) and (-Delta)^T = Delta.
G_T = expm(Delta)

# Step 3: Compute the post-reset covariance matrix, Sigma_prime.
Sigma_prime = G @ Sigma @ G_T

# Step 4: Print all components of the final equation.
print("--- Attitude Reset Covariance Calculation ---")
print("\nInitial attitude deviation vector (delta):\n", delta_vec)
print("\nPre-reset covariance matrix (Sigma):\n", Sigma)
print("\n---------------------------------------------")
print("\nThe final equation is: Sigma' = exp(-Delta) * Sigma * exp(Delta)")
print("\nWhere:")

print("\nexp(-Delta) = \n")
print(G)
print("\nSigma = \n")
print(Sigma)
print("\nexp(Delta) = \n")
print(G_T)

print("\nResulting post-reset covariance matrix (Sigma'):\n")
print(Sigma_prime)