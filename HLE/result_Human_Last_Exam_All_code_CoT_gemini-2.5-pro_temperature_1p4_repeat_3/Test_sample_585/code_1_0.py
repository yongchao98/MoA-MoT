import numpy as np

def hat(v):
    """
    Computes the skew-symmetric matrix for a 3-element vector v.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def rodrigues(phi):
    """
    Computes the rotation matrix R = exp(hat(phi)) using Rodrigues' formula.
    """
    phi = np.asarray(phi)
    angle = np.linalg.norm(phi)
    if angle < 1e-12:
        return np.identity(3)
    
    axis = phi / angle
    K = hat(axis)
    # Rodrigues' rotation formula
    R = np.identity(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    return R

def right_jacobian_so3(phi):
    """
    Computes the right Jacobian of SO(3).
    """
    phi = np.asarray(phi)
    angle = np.linalg.norm(phi)
    K = hat(phi)
    
    if angle < 1e-9:
        # Use Taylor series expansion for small angles
        return np.identity(3) - 0.5 * K + (1./6.) * (K @ K)
    
    K_sq = K @ K
    a = (1 - np.cos(angle)) / (angle**2)
    b = (angle - np.sin(angle)) / (angle**3)
    
    Jr = np.identity(3) - a * K + b * K_sq
    return Jr

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance matrix.
    
    Args:
        delta (np.ndarray): The 3x1 attitude deviation vector.
        Sigma (np.ndarray): The 3x3 pre-reset covariance matrix.
        
    Returns:
        np.ndarray: The 3x3 post-reset covariance matrix.
    """
    # Compute the rotation matrix R = exp(hat(delta))
    R_delta = rodrigues(delta)
    
    # Compute the right Jacobian of SO(3)
    Jr_delta = right_jacobian_so3(delta)
    
    # The transformation Jacobian J for the covariance update
    J = R_delta @ Jr_delta
    
    # Compute the post-reset covariance
    Sigma_prime = J @ Sigma @ J.T
    
    return Sigma_prime

if __name__ == '__main__':
    # --- Example Usage ---
    # Estimated attitude deviation vector (radians)
    # This is the 'delta' used to update the reference attitude
    delta = np.array([0.1, -0.2, 0.05])

    # Pre-reset covariance matrix
    # This represents the uncertainty of the attitude deviation estimate
    Sigma = np.diag([
        np.deg2rad(0.5)**2,
        np.deg2rad(0.5)**2,
        np.deg2rad(1.0)**2
    ])
    
    # Compute the post-reset covariance
    Sigma_prime = compute_post_reset_covariance(delta, Sigma)
    
    # The final expression is Sigma' = J * Sigma * J^T
    # where J = exp(hat(delta)) * Jr(delta)
    R_delta_calc = rodrigues(delta)
    Jr_delta_calc = right_jacobian_so3(delta)
    J_calc = R_delta_calc @ Jr_delta_calc

    print("--- Inputs ---")
    print(f"Attitude deviation delta = {delta}")
    print("\nPre-reset covariance Sigma = \n", Sigma)
    
    print("\n--- Calculation ---")
    print("The post-reset covariance Sigma' is calculated using the formula:")
    print("Sigma' = (exp(hat(delta)) * Jr(delta)) * Sigma * (exp(hat(delta)) * Jr(delta))^T")
    print("\nWhere:")
    print("exp(hat(delta)) = \n", R_delta_calc)
    print("\nJr(delta) = \n", Jr_delta_calc)
    print("\nTransformation Jacobian J = exp(hat(delta)) * Jr(delta) = \n", J_calc)

    print("\n--- Final Result ---")
    print("Post-reset covariance Sigma' = \n", Sigma_prime)

    final_expression_str = f"Sigma' = \n{J_calc} \n* \n{Sigma} \n* \n{J_calc.T}"
    
    print("\nThe full equation with the final numbers is:")
    print(final_expression_str)
