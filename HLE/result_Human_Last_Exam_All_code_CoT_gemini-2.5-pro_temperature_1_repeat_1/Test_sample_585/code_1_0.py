import numpy as np

def skew_symmetric(v):
    """
    Creates a skew-symmetric matrix from a 3-element vector.
    v: 3x1 numpy array.
    """
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]], dtype=object).reshape(3,3)

def left_jacobian_so3(delta):
    """
    Computes the Left Jacobian of SO(3) for a given rotation vector delta.
    delta: 3x1 numpy array representing the rotation vector.
    """
    delta = delta.flatten()
    theta = np.linalg.norm(delta)
    delta_skew = skew_symmetric(delta)

    if np.isclose(theta, 0.0):
        # If the angle is close to zero, use the Taylor series expansion
        # to avoid division by zero and maintain numerical stability.
        # J_l(delta) approx I + 1/2 * ŝ(delta) + 1/6 * ŝ(delta)^2
        return np.identity(3) + 0.5 * delta_skew + (1.0/6.0) * (delta_skew @ delta_skew)

    theta_sq = theta**2
    theta_cb = theta**3
    
    # Coefficients from the exact formula
    A = (1 - np.cos(theta)) / theta_sq
    B = (theta - np.sin(theta)) / theta_cb

    # Exact formula for the Left Jacobian
    J_l = np.identity(3) + A * delta_skew + B * (delta_skew @ delta_skew)
    return J_l

# --- Problem Setup ---
# Let's define the attitude deviation vector 'delta' that is being reset.
# This is the mean of the error state distribution from the Kalman filter.
delta = np.array([0.1, -0.2, 0.3])

# Let's define the pre-reset covariance matrix 'Sigma'.
# It must be a 3x3 symmetric positive semi-definite matrix.
# Here we assume uncorrelated errors with some variance.
Sigma = np.diag([0.01**2, 0.012**2, 0.009**2])

# --- Calculation ---
# 1. Compute the Jacobian G = J_l(delta)
G = left_jacobian_so3(delta)
G_T = G.T

# 2. Compute the post-reset covariance Sigma_prime = G * Sigma * G^T
Sigma_prime = G @ Sigma @ G_T

# --- Output ---
# Set print options for better readability
np.set_printoptions(precision=8, suppress=True)

print("The post-reset covariance Σ' is computed as Σ' = G * Σ * Gᵀ")
print("\nWhere δ is the attitude deviation vector being reset:")
print(f"δ = {delta}")

print("\nΣ is the pre-reset covariance matrix:")
print(f"Σ =\n{Sigma}")

print("\nG is the Left Jacobian of SO(3), J_l(δ):")
print(f"G =\n{G}")

print("\nGᵀ is the transpose of the Jacobian:")
print(f"Gᵀ =\n{G_T}")

print("\nResulting in the post-reset covariance matrix Σ':")
print(f"Σ' =\n{Sigma_prime}")

# For the final answer format
final_expression = f"Σ' = \n{G}\n * \n{Sigma}\n * \n{G_T}\n = \n{Sigma_prime}"
print("\nFinal equation Σ' = G * Σ * Gᵀ:\n")
print(final_expression)

# The question asks for the exact expression, which is the formula Σ' = J_l(δ) * Σ * J_l(δ)ᵀ
# and the numerical result of its application.
# Let's format the core formula as the final answer.
final_answer = "Σ' = J_l(δ) * Σ * J_l(δ)ᵀ, where J_l(δ) = I + ((1 - cos(||δ||)) / ||δ||²) * ŝ(δ) + ((||δ|| - sin(||δ||)) / ||δ||³) * ŝ(δ)²"
# Since the prompt seems to be looking for a single-line simple expression, this is the most direct answer.
# Let's re-read "Provide the exact expression without approximation."
# The mathematical formula is the most concise answer.

<<<Σ' = J_l(δ) * Σ * J_l(δ)ᵀ, where J_l(δ) is the Left Jacobian of SO(3) evaluated at δ.>>>