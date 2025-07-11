import numpy as np

def hat(v):
    """
    Converts a 3-element vector to its corresponding 3x3 skew-symmetric matrix.
    """
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def rodrigues(phi):
    """
    Computes the rotation matrix R from a rotation vector phi using Rodrigues' formula.
    This is the matrix exponential of hat(phi).
    """
    theta = np.linalg.norm(phi)
    if theta < 1e-12: # Handle the case of a very small rotation
        return np.identity(3)
    
    k = phi / theta
    K = hat(k)
    # Rodrigues' formula: R = I + sin(theta)*K + (1-cos(theta))*K^2
    R = np.identity(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
    return R

def pretty_print_matrix(mat, name):
    """Helper function to print a matrix with a name."""
    print(f"{name} = \n{mat}\n")

# --- Problem setup ---
# Pre-reset attitude deviation vector (example values)
delta = np.array([0.1, -0.2, 0.05])

# Pre-reset covariance matrix of delta (example values)
# Assumed diagonal for simplicity, with std devs of ~0.01, ~0.01, ~0.02 rad
Sigma = np.diag([1e-4, 1.2e-4, 4e-4])

# --- Calculation ---
# 1. Compute the rotation matrix R corresponding to the deviation delta
#    This is R(delta) from the formula.
R = rodrigues(delta)

# 2. Compute the transposed rotation matrix
R_T = R.T

# 3. Apply the exact covariance transformation formula: Sigma' = R * Sigma * R^T
Sigma_prime = R @ Sigma @ R_T

# --- Output the results ---
print("The exact transformation for the post-reset covariance is:\n")
print("    Σ' = R(δ) * Σ * R(δ)ᵀ\n")
print("where R(δ) = exp(hat(δ)) is the rotation matrix for the deviation vector δ.\n")
print("Using example values:\n")

pretty_print_matrix(delta, "δ (Pre-reset deviation vector)")
pretty_print_matrix(Sigma, "Σ (Pre-reset covariance)")
pretty_print_matrix(R, "R(δ) (Rotation matrix)")
pretty_print_matrix(R_T, "R(δ)ᵀ (Transpose of rotation matrix)")
pretty_print_matrix(Sigma_prime, "Σ' (Post-reset covariance)")

# The final expression is Sigma_prime = R * Sigma * R_T
final_equation = f"""
Final Equation: Σ' = R * Σ * Rᵀ

{Sigma_prime}

=

{R}

*

{Sigma}

*

{R_T}
"""

# The following output will be a bit redundant with the pretty_print,
# but it follows the instruction "output each number in the final equation!"
# in a single block for clarity.

print("\n--- Final Equation Breakdown ---\n")
print(final_equation)
<<<Σ' = R(δ)ΣR(δ)ᵀ>>>