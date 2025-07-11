import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix form.
    hat(v) = [[ 0, -v3,  v2],
              [ v3,   0, -v1],
              [-v2,  v1,   0]]
    """
    if not isinstance(v, np.ndarray) or v.shape != (3,):
        raise ValueError("Input must be a 3-element NumPy vector.")
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def rodrigues(v):
    """
    Computes the matrix exponential exp(hat(v)) of a 3D vector v using Rodrigues' formula.
    This gives the rotation matrix corresponding to the rotation vector v.
    """
    if not isinstance(v, np.ndarray) or v.shape != (3,):
        raise ValueError("Input must be a 3-element NumPy vector.")

    theta = np.linalg.norm(v)
    
    if theta < 1e-12: # Handle the case of a very small rotation
        return np.identity(3)
        
    k = v / theta
    K = skew(k)
    
    # Rodrigues' formula
    # R = I + sin(theta)*K + (1-cos(theta))*K^2
    R = np.identity(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
    return R

def pretty_print(name, matrix):
    """Helper function for clean printing of matrices."""
    print(f"{name}:\n{matrix}\n")

# --- Problem Setup ---
# Let's assume the Kalman filter has estimated a non-zero attitude deviation `delta`
# with an associated covariance `Sigma`.

# 1. The attitude deviation vector delta (δ)
delta = np.array([0.1, -0.2, 0.05])

# 2. The covariance matrix Sigma (Σ) associated with delta
# (must be symmetric and positive semi-definite)
Sigma = np.array([
    [0.01, 0.001, -0.002],
    [0.001, 0.02, 0.005],
    [-0.002, 0.005, 0.015]
])

print("--- Attitude Reset Covariance Calculation ---\n")
pretty_print("Initial deviation vector δ", delta)
pretty_print("Initial covariance matrix Σ", Sigma)
print("The post-reset covariance Σ' is calculated using the formula:")
print("Σ' = exp(hat(-δ)) * Σ * exp(hat(δ))\n")
print("--- Step-by-step Calculation ---\n")

# --- Calculation ---

# 1. Calculate the rotation matrix from the deviation, R_delta = exp(hat(δ))
R_delta = rodrigues(delta)
pretty_print("exp(hat(δ))", R_delta)

# 2. Calculate the rotation matrix from the negative deviation, exp(hat(-δ))
# This is also the transpose of R_delta.
R_minus_delta = rodrigues(-delta)
pretty_print("exp(hat(-δ))", R_minus_delta)

# 3. Compute the new covariance matrix Sigma_prime (Σ')
# Σ' = exp(hat(-δ)) * Σ * exp(hat(δ))
Sigma_prime = R_minus_delta @ Sigma @ R_delta

print("--- Final Result ---\n")
pretty_print("Post-reset covariance matrix Σ'", Sigma_prime)

# The final expression is Sigma_prime = exp(hat(-delta)) * Sigma * exp(hat(delta))
# where each term is printed above.
