import numpy as np

def skew(v):
    """
    Converts a 3-element vector to a skew-symmetric matrix.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def calculate_post_reset_covariance(delta_hat, Sigma):
    """
    Calculates the post-reset covariance matrix for an attitude deviation.

    Args:
        delta_hat (np.ndarray): The 3-element attitude deviation vector estimate (mean).
        Sigma (np.ndarray): The 3x3 pre-reset covariance matrix of the deviation.

    Returns:
        np.ndarray: The 3x3 post-reset covariance matrix.
    """
    # The transformation Jacobian is the left Jacobian of SO(3), J_l(-delta_hat).
    phi = -delta_hat
    angle = np.linalg.norm(phi)
    
    phi_skew = skew(phi)
    phi_skew_sq = phi_skew @ phi_skew
    
    # Use Taylor series expansion for small angles for numerical stability
    if angle < 1e-8:
        A = 1/2
        B = 1/6
    else:
        angle_sq = angle**2
        angle_cub = angle**3
        A = (1 - np.cos(angle)) / angle_sq
        B = (angle - np.sin(angle)) / angle_cub
        
    # The Jacobian G = J_l(-delta_hat)
    G = np.eye(3) + A * phi_skew + B * phi_skew_sq
    
    # The post-reset covariance is Sigma' = G * Sigma * G^T
    Sigma_prime = G @ Sigma @ G.T
    
    return G, Sigma_prime

if __name__ == '__main__':
    # Example pre-reset values
    # The estimated attitude deviation (mean)
    delta_hat = np.array([0.1, -0.2, 0.15])
    
    # The pre-reset covariance matrix
    Sigma = np.array([
        [0.01, 0.001, -0.002],
        [0.001, 0.012, 0.003],
        [-0.002, 0.003, 0.015]
    ])
    
    # Calculate the post-reset covariance
    G, Sigma_prime = calculate_post_reset_covariance(delta_hat, Sigma)
    
    print("This script calculates the post-reset covariance matrix Σ' for an attitude deviation.")
    print("\nThe formula is: Σ' = G * Σ * G^T, where G = J_l(-δ̂)\n")

    print("Given δ̂ (delta_hat):")
    print(delta_hat)
    print("\nAnd pre-reset covariance Σ (Sigma):")
    print(Sigma)
    
    print("\nThe computed Jacobian G = J_l(-δ̂) is:")
    print(G)
    
    # The final equation's result
    print("\nThe final post-reset covariance matrix Σ' is:")
    print(Sigma_prime)
    
    final_equation_string = f"\nΣ' =\n{np.array_str(Sigma_prime, precision=5)}"
    
    print("\nEach number in the final equation:")
    print(final_equation_string)

<<<Σ' = J_l(-δ̂) * Σ * J_l(-δ̂)ᵀ>>>