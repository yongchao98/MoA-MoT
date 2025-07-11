import numpy as np

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance matrix for an attitude deviation.

    Args:
        delta (np.ndarray): The 3x1 attitude deviation vector used for the reset.
        Sigma (np.ndarray): The 3x3 pre-reset covariance matrix of the deviation.

    Returns:
        np.ndarray: The 3x3 post-reset covariance matrix Sigma_prime.
    """
    # Ensure delta is a 3-element array
    delta = np.asarray(delta).flatten()
    if delta.shape[0] != 3:
        raise ValueError("delta must be a 3-element vector.")
    
    # Calculate the norm of the deviation vector
    theta = np.linalg.norm(delta)
    
    # Construct the skew-symmetric matrix of delta
    delta_x = np.array([
        [0, -delta[2], delta[1]],
        [delta[2], 0, -delta[0]],
        [-delta[1], delta[0], 0]
    ])
    
    # Compute the Jacobian matrix G = J_r(delta)
    # To avoid division by zero for small theta, we use Taylor series expansions
    if theta < 1e-8:
        # For theta -> 0, G -> I (Identity matrix)
        # Using Taylor expansion: A -> 1/2, B -> 1/6
        A = 1/2
        B = 1/6
    else:
        theta2 = theta**2
        theta3 = theta**3
        A = (1 - np.cos(theta)) / theta2
        B = (theta - np.sin(theta)) / theta3
        
    I = np.eye(3)
    delta_x2 = delta_x @ delta_x
    G = I + A * delta_x + B * delta_x2  # Note: some literature uses G = I - A*delta_x + B*delta_x2, which is Jr(-delta). The convention here corresponds to δ_new = δ_old - δ_est
    
    # Correct Jacobian should be J_r(-delta) for the transformation δ' = δ - δ_est
    # J_r(-delta) = I - A * [δ]x + B * [δ]x^2
    # Many texts establish that G = I - [δ]x/2 for small δ. To get the sign right:
    # A for small theta is 1/2. G = I + 1/2 [δ]x + ...
    # Wait, the relationship exp(δ_true) = exp(δ_est) exp(δ_new) leads to 
    # the jacobian being G = J_r(-δ_est). Let's use that one, it's more standard.
    # J_r(-x) = I - (1-cos)/||x||^2 [x]x + (||x||-sin)/||x||^3 [x]x^2
    G = I - A * delta_x + B * delta_x2

    # Propagate the covariance: Sigma_prime = G * Sigma * G^T
    Sigma_prime = G @ Sigma @ G.T
    
    # --- Output the equation as requested ---
    print("The post-reset covariance Σ' is computed as: Σ' = G * Σ * G^T\n")
    
    print("The pre-reset attitude deviation vector δ was:")
    print(delta)
    print("\nThe pre-reset covariance matrix Σ was:")
    print(Sigma)
    
    print("\nResulting Jacobian G = J_r(-δ):")
    print(G)
    
    print("\nTranspose of the Jacobian, G^T:")
    print(G.T)
    
    print("\nFinal post-reset covariance matrix Σ' is:")
    print(Sigma_prime)

    return Sigma_prime

if __name__ == '__main__':
    # Define an example deviation vector (the estimate used for the reset)
    # Represents a rotation of approx. 0.26 radians around axis ~[0.4, 0.8, -0.5]
    delta_vector = np.array([0.1, 0.2, -0.15])
    
    # Define an example pre-reset covariance matrix
    # (e.g., indicating uncertainty of 0.1 rad^2 on each axis)
    pre_reset_sigma = np.diag([0.01, 0.01, 0.01])
    
    np.set_printoptions(precision=6, suppress=True)
    
    compute_post_reset_covariance(delta_vector, pre_reset_sigma)