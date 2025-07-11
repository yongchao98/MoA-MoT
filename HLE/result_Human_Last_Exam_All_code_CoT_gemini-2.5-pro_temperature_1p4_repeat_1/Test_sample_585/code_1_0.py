import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix form.
    """
    if not hasattr(v, '__len__') or len(v) != 3:
        raise ValueError("Input must be a 3-element vector.")
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def compute_post_reset_covariance(delta_hat, Sigma):
    """
    Computes the post-reset covariance matrix Σ' using the exact formula.

    Args:
        delta_hat (np.ndarray): The 3x1 attitude deviation vector used for the reset.
        Sigma (np.ndarray): The 3x3 covariance matrix of the deviation before the reset.

    Returns:
        A tuple containing (J_r, Sigma_prime)
    """
    # Ensure inputs are numpy arrays
    delta_hat = np.asarray(delta_hat).flatten()
    Sigma = np.asarray(Sigma)

    # Calculate the magnitude of the rotation vector
    theta = np.linalg.norm(delta_hat)
    
    # Skew-symmetric matrix of delta_hat
    S = skew(delta_hat)
    
    # Identity matrix
    I = np.identity(3)
    
    # Use Taylor series expansion for small angles to avoid division by zero
    if theta < 1e-8:
        # A -> -1/2, B -> 1/6 for small theta
        A = -0.5
        B = 1.0 / 6.0
        # For very small angles, Jr can be approximated more simply
        J_r = I + A * S + B @ (S @ S)

    else:
        theta2 = theta * theta
        theta3 = theta2 * theta
        
        # Coefficients for the Jacobian formula
        # The standard formula for Jr uses a '-' before the first term, so our coefficient A is positive.
        A = (1 - np.cos(theta)) / theta2
        B = (theta - np.sin(theta)) / theta3

        # Exact formula for the Right Jacobian of SO(3)
        J_r = I - A * S + B * (S @ S)

    # Propagate the covariance
    # Sigma_prime = J_r * Sigma * J_r^T
    Sigma_prime = J_r @ Sigma @ J_r.T
    
    return J_r, Sigma_prime

if __name__ == '__main__':
    # --- Example Inputs ---
    # The attitude deviation vector used in the reset step (e.g., in radians)
    delta_hat = np.array([0.1, -0.05, 0.08])
    
    # The covariance of the attitude deviation before the reset
    # (e.g., corresponding to 0.5 degrees standard deviation on each axis)
    std_dev_rad = np.deg2rad(0.5)
    Sigma = np.diag([std_dev_rad**2, std_dev_rad**2, std_dev_rad**2])

    # --- Calculation ---
    J_r, Sigma_prime = compute_post_reset_covariance(delta_hat, Sigma)

    # --- Output Results ---
    np.set_printoptions(precision=8, suppress=True)

    print("The post-reset covariance Σ' is computed as: Σ' = J * Σ * Jᵀ\n")
    
    print("Pre-reset deviation δ:")
    print(delta_hat)
    print("\nPre-reset covariance Σ:")
    print(Sigma)

    print("\nComputed Right Jacobian J = J_r(δ):")
    print(J_r)

    print("\nJacobian Transposed Jᵀ:")
    print(J_r.T)
    
    print("\nThe full equation is:")
    print("Σ' = \n", J_r, "\n * \n", Sigma, "\n * \n", J_r.T)
    
    print("\nResulting Post-reset covariance Σ':")
    print(Sigma_prime)

    # The final expression itself is the formula applied in the code
    # Final answer format for the symbolic expression
    final_answer = "Σ' = J_r(δ) * Σ * J_r(δ)ᵀ, where J_r(δ) = I - ((1-cos(θ))/θ²) * skew(δ) + ((θ-sin(θ))/θ³) * skew(δ)² and θ=||δ||"
    # For a numerical answer, we provide the final matrix from the calculation
    # We will output this symbolic expression since it is the core of the answer
    
    print(f"\n\n<<<Σ' = J_r(δ) * Σ * J_r(δ)ᵀ>>>")
