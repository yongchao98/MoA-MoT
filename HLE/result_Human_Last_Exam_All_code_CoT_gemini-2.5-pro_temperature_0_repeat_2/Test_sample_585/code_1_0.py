import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix form.
    """
    if not isinstance(v, np.ndarray) or v.shape != (3,):
        raise ValueError("Input must be a 3-element numpy array.")
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def left_jacobian_so3(phi):
    """
    Computes the Left Jacobian of SO(3) for a given rotation vector phi.
    This is the exact expression without approximation.
    """
    if not isinstance(phi, np.ndarray) or phi.shape != (3,):
        raise ValueError("Input must be a 3-element numpy array.")
        
    angle = np.linalg.norm(phi)
    
    # Handle the case where the angle is close to zero to avoid division by zero.
    # For very small angles, the Jacobian is close to the identity matrix.
    if np.isclose(angle, 0.0):
        return np.identity(3)
        
    axis = phi / angle
    K = skew(axis)
    K_sq = K @ K
    
    # Exact formula for the Left Jacobian of SO(3)
    # J_l(phi) = I + (1 - cos(angle)) / angle**2 * hat(phi) + (angle - sin(angle)) / angle**3 * hat(phi)**2
    # Since hat(phi) = angle * K, we can substitute to get:
    # J_l(phi) = I + (1 - cos(angle)) / angle * K + (angle - sin(angle)) / angle**2 * K_sq
    # This form is more numerically stable.
    
    c1 = (1 - np.cos(angle)) / angle
    c2 = (angle - np.sin(angle)) / (angle**2)
    
    J = np.identity(3) + c1 * skew(phi) + c2 * (skew(phi) @ skew(phi))
    
    return J

def main():
    """
    Main function to demonstrate the covariance reset calculation.
    """
    # --- Inputs ---
    # A sample attitude deviation vector 'delta' (in radians)
    delta = np.array([0.1, -0.2, 0.15])
    
    # A sample pre-reset covariance matrix 'Sigma' for the deviation vector.
    # (e.g., corresponding to standard deviations of ~0.5 degrees)
    std_dev = np.deg2rad(0.5)
    Sigma = np.diag([std_dev**2, std_dev**2, std_dev**2])

    # --- Calculation ---
    # 1. Compute the Left Jacobian matrix J_l(delta)
    J = left_jacobian_so3(delta)
    
    # 2. Compute the post-reset covariance matrix Sigma'
    # Sigma' = J * Sigma * J^T
    Sigma_prime = J @ Sigma @ J.T
    
    # --- Output ---
    np.set_printoptions(precision=8, suppress=True)
    
    print("The post-reset covariance Σ' is computed as: Σ' = J * Σ * Jᵀ\n")
    
    print("Given the pre-reset deviation δ:")
    print(delta)
    print("\nAnd the pre-reset covariance Σ:")
    print(Sigma)
    
    print("\nFirst, we compute the transformation matrix J = J_l(δ):")
    print(J)
    
    print("\nThen, we compute its transpose Jᵀ:")
    print(J.T)
    
    print("\nFinally, the post-reset covariance Σ' = J * Σ * Jᵀ is:")
    print(Sigma_prime)

if __name__ == "__main__":
    main()