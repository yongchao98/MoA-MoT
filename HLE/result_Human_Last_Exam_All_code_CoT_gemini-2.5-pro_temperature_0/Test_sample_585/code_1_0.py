import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix.
    """
    v = v.flatten()
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance matrix for an attitude deviation.

    Args:
        delta (np.ndarray): The 3x1 attitude deviation vector used for the reset.
        Sigma (np.ndarray): The 3x3 covariance matrix of the attitude deviation.

    Returns:
        np.ndarray: The 3x3 post-reset covariance matrix Sigma_prime.
    """
    print("--- Input Parameters ---")
    print(f"δ (delta vector):\n{delta.flatten()}")
    print(f"Σ (pre-reset covariance):\n{Sigma}\n")

    phi = np.linalg.norm(delta)
    
    # Handle the case of a very small rotation to avoid division by zero
    if phi < 1e-9:
        print("Delta is near zero. Post-reset covariance is the same as pre-reset covariance.")
        return Sigma

    delta_hat = skew(delta)
    delta_hat_sq = delta_hat @ delta_hat
    
    # Coefficients for R and J_l
    phi2 = phi * phi
    phi3 = phi2 * phi
    
    s_phi = np.sin(phi)
    c_phi = np.cos(phi)
    
    # Coefficients for Rotation Matrix R
    c1_R = s_phi / phi
    c2_R = (1 - c_phi) / phi2
    
    # Coefficients for Left Jacobian J_l
    c2_Jl = (1 - c_phi) / phi2
    c3_Jl = (phi - s_phi) / phi3
    
    I = np.identity(3)
    
    # Compute R(δ) and J_l(δ)
    R = I + c1_R * delta_hat + c2_R * delta_hat_sq
    J_l = I + c2_Jl * delta_hat + c3_Jl * delta_hat_sq
    
    # The transformation matrix G = R(δ) * J_l(δ)
    G = R @ J_l
    
    # Compute the post-reset covariance Σ' = G * Σ * G^T
    Sigma_prime = G @ Sigma @ G.T
    
    print("--- Intermediate Computations ---")
    print(f"φ (magnitude of delta): {phi:.4f}")
    print(f"R(δ) (Rotation Matrix):\n{R}\n")
    print(f"J_l(δ) (Left Jacobian):\n{J_l}\n")
    print(f"G = R(δ) * J_l(δ) (Transformation Matrix):\n{G}\n")
    
    print("--- Final Equation ---")
    print("Σ' = G * Σ * G^T")
    print("Where G, Σ are the matrices above.\n")
    
    print("--- Result ---")
    print(f"Σ' (post-reset covariance):\n{Sigma_prime}")
    
    return Sigma_prime

if __name__ == '__main__':
    # Example usage:
    # Define an attitude deviation vector δ (e.g., in radians)
    delta_k = np.array([0.1, 0.2, -0.1])

    # Define the covariance matrix Σ associated with δ
    # (e.g., uncertainty of 0.01 rad on each axis)
    sigma_rad = 0.01
    Sigma_k = np.diag([sigma_rad**2, sigma_rad**2, sigma_rad**2])

    # Compute the post-reset covariance
    Sigma_k_prime = compute_post_reset_covariance(delta_k, Sigma_k)
