import numpy as np

def skew(v):
    """
    Converts a 3-element vector to a skew-symmetric matrix.
    v_hat = skew(v) such that v_hat @ w = v x w.
    """
    if not isinstance(v, np.ndarray) or v.shape != (3,):
        raise ValueError("Input must be a 3-element numpy array.")
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def rodrigues_exp(v):
    """
    Computes the matrix exponential exp(hat(v)) for a 3D vector v (rotation vector)
    using Rodrigues' formula.
    """
    theta = np.linalg.norm(v)
    if theta < 1e-12:  # Treat as zero for numerical stability
        return np.identity(3)
    
    v_hat = skew(v)
    v_hat_sq = v_hat @ v_hat
    
    c1 = np.sin(theta) / theta
    c2 = (1 - np.cos(theta)) / (theta**2)
    
    R = np.identity(3) + c1 * v_hat + c2 * v_hat_sq
    return R

def left_jacobian(v):
    """
    Computes the left Jacobian of SO(3), J_l(v).
    """
    theta = np.linalg.norm(v)
    if theta < 1e-12: # Treat as zero for numerical stability
        return np.identity(3)
        
    v_hat = skew(v)
    v_hat_sq = v_hat @ v_hat
    
    c1 = (1 - np.cos(theta)) / (theta**2)
    c2 = (theta - np.sin(theta)) / (theta**3)
    
    J = np.identity(3) + c1 * v_hat + c2 * v_hat_sq
    return J

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance Σ' using the exact, non-approximated formula.
    
    Args:
        delta (np.array): The 3x1 attitude deviation vector being reset.
        Sigma (np.array): The 3x3 covariance matrix associated with delta.
    
    Returns:
        np.array: The 3x3 post-reset covariance matrix Σ'.
    """
    print("--- Input Values ---")
    print(f"δ (delta vector):\n{delta}\n")
    print(f"Σ (pre-reset covariance):\n{Sigma}\n")

    # --- Calculation Steps ---
    
    # 1. Rotation matrix from delta
    R = rodrigues_exp(delta)
    
    # 2. Left Jacobian of SO(3)
    Jl = left_jacobian(delta)
    
    # 3. The full transformation Jacobian G
    G = R @ Jl
    
    # 4. The post-reset covariance Sigma' = G * Sigma * G^T
    Sigma_prime = G @ Sigma @ G.T
    
    print("--- Intermediate Matrices for the Equation: Σ' = G * Σ * G^T ---")
    print(f"exp(hat(δ)) (Rotation Matrix R):\n{R}\n")
    print(f"J_l(δ) (Left Jacobian):\n{Jl}\n")
    print(f"G = exp(hat(δ)) * J_l(δ) (Full Transformation Jacobian):\n{G}\n")

    print("--- Final Result ---")
    print(f"Σ' (post-reset covariance):\n{Sigma_prime}\n")

    return Sigma_prime

if __name__ == '__main__':
    # Define an example attitude deviation vector `delta` and its covariance `Sigma`
    # `delta` is the estimated error that will be absorbed into the reference attitude
    delta_vector = np.array([0.1, -0.15, 0.2])

    # `Sigma` is the covariance of the true error `delta_true`
    Sigma_matrix = np.array([
        [0.01, 0.001, -0.002],
        [0.001, 0.012, 0.003],
        [-0.002, 0.003, 0.015]
    ])
    
    print("The exact expression for the post-reset covariance Σ' is:")
    print("Σ' = exp(hat(δ)) * J_l(δ) * Σ * J_l(δ)^T * exp(hat(δ))^T\n")
    print("where δ is the reset vector, Σ is the pre-reset covariance,")
    print("exp(hat(δ)) is the matrix exponential, and J_l(δ) is the Left Jacobian of SO(3).\n")
    print("Here is a numerical example:\n")

    # Compute and print the results
    post_reset_cov = compute_post_reset_covariance(delta_vector, Sigma_matrix)

    # Final answer format for the value of Sigma_prime
    # We will format the final numpy array as a string to fit the required format.
    final_answer_string = np.array2string(post_reset_cov, separator=', ', formatter={'float_kind':lambda x: f"{x:.6f}"})
    # Remove newline characters and extra spaces for a compact representation
    final_answer_string = final_answer_string.replace('\n', '').replace(' ', '')
    
    print("<<<" + final_answer_string + ">>>")
