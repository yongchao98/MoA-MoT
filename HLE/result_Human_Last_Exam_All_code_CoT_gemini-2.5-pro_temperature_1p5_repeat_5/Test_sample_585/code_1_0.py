import numpy as np

def skew(v):
    """
    Creates a skew-symmetric matrix from a 3-element vector.
    v_skew = [[  0, -v2,  v1],
              [ v2,   0, -v0],
              [-v1,  v0,   0]]
    """
    v = np.asarray(v)
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def right_jacobian_so3(phi):
    """
    Computes the right Jacobian of SO(3).

    Jr(phi) = I - (1-cos(||phi||))/||phi||^2 * phi_hat + 
              (||phi|| - sin(||phi||))/||phi||^3 * phi_hat^2
    """
    phi = np.asarray(phi)
    norm_phi = np.linalg.norm(phi)
    
    # Handle the case where the angle is very small to avoid division by zero
    if norm_phi < 1e-9:
        # For phi -> 0, Jr(phi) -> I. We can use a Taylor expansion for more accuracy if needed.
        # Jr(phi) approx I - 1/2*phi_hat + 1/6*phi_hat^2
        return np.identity(3)

    phi_hat = skew(phi)
    phi_hat_sq = phi_hat @ phi_hat
    
    A = (1 - np.cos(norm_phi)) / (norm_phi**2)
    B = (norm_phi - np.sin(norm_phi)) / (norm_phi**3)
    
    # The formula for Jr(phi)
    J = np.identity(3) - A * phi_hat + B * phi_hat_sq
    return J

def calculate_post_reset_covariance(delta, Sigma):
    """
    Calculates the post-reset covariance Sigma' = J * Sigma * J^T,
    where J = Jr(-delta).
    """
    # The Jacobian of the reset transformation is the right Jacobian of SO(3)
    # evaluated at -delta.
    J = right_jacobian_so3(-delta)

    # The post-reset covariance using the EKF propagation rule
    Sigma_prime = J @ Sigma @ J.T

    # --- Output the results ---
    print("The exact expression for the post-reset covariance is: Σ' = J * Σ * Jᵀ")
    print("where J is the exact Jacobian of the attitude error transformation, J = Jᵣ(-δ).\n")
    
    np.set_printoptions(precision=6, suppress=True)

    print(f"Given reset vector δ =\n{delta}\n")
    print(f"And pre-reset covariance Σ =\n{Sigma}\n")

    print(f"The Jacobian J = Jᵣ(-δ) is calculated as:\n{J}\n")

    print("--- Final Equation Components ---")
    print("Σ' = \n")
    print(J)
    print("\n * \n")
    print(Sigma)
    print("\n * \n")
    print(J.T)
    print("\n = \n")
    print(Sigma_prime)
    
    return Sigma_prime

if __name__ == '__main__':
    # Define a sample reset vector δ (delta)
    # This would be the filter's estimate of the attitude deviation
    delta_vector = np.array([0.1, -0.05, 0.2])

    # Define a sample pre-reset covariance matrix Σ (Sigma)
    # This is the covariance associated with the deviation state δ
    Sigma_matrix = np.diag([0.01**2, 0.01**2, 0.02**2])

    # Calculate and print the post-reset covariance
    post_reset_Sigma = calculate_post_reset_covariance(delta_vector, Sigma_matrix)

    # For verification, the final result is also returned
    final_answer = post_reset_Sigma
    # The output format is not suitable for a matrix, so we will not use the <<<>>> format.
