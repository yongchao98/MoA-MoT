import numpy as np

def skew(v):
    """
    Converts a 3-element vector to a 3x3 skew-symmetric matrix (hat operator).
    
    Args:
        v (np.ndarray): A 3-element vector.
        
    Returns:
        np.ndarray: The corresponding 3x3 skew-symmetric matrix.
    """
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def exp_map(v):
    """
    Computes the rotation matrix from a rotation vector (exponential map for SO(3))
    using Rodrigues' formula.
    
    Args:
        v (np.ndarray): A 3-element rotation vector.
        
    Returns:
        np.ndarray: The corresponding 3x3 rotation matrix.
    """
    theta = np.linalg.norm(v)
    K = skew(v)
    if theta < 1e-9:  # Use Taylor series expansion for small angles to avoid division by zero
        # R = I + K + 0.5*K^2
        return np.identity(3) + K + 0.5 * (K @ K)
    else:
        # R = I + (sin(theta)/theta)*K + ((1-cos(theta))/theta^2)*K^2
        return np.identity(3) + (np.sin(theta) / theta) * K + ((1 - np.cos(theta)) / (theta**2)) * (K @ K)

def inv_left_jacobian(v):
    """
    Computes the inverse of the left Jacobian of SO(3).
    
    Args:
        v (np.ndarray): A 3-element rotation vector.
        
    Returns:
        np.ndarray: The 3x3 inverse left Jacobian matrix.
    """
    theta = np.linalg.norm(v)
    K = skew(v)
    if theta < 1e-9: # Use Taylor series expansion for small angles
        # J_l(v)^-1 = I + 0.5*K + (1/12)*K^2
        return np.identity(3) + 0.5 * K + (1.0/12.0) * (K @ K)
    else:
        # Full formula using cotangent for better numerical stability near multiples of pi
        # J_l(v)^-1 = I + 0.5*K + (1/theta^2 - cot(theta/2)/(2*theta)) * K^2
        theta2 = theta**2
        half_theta = theta / 2.0
        cot_half_theta = 1.0 / np.tan(half_theta)
        # Coefficient for the K^2 term
        coef = (1.0 / theta2) - (cot_half_theta / (2.0 * theta))
        return np.identity(3) + 0.5 * K + coef * (K @ K)

def main():
    """
    Main function to demonstrate the covariance reset calculation.
    """
    # The attitude deviation vector `delta` that is being reset.
    # This represents the mean of the attitude error distribution before the reset.
    delta = np.array([0.1, -0.2, 0.3])

    # The pre-reset covariance matrix `Sigma` associated with the error `e = delta_true - delta`.
    # We assume an initial uncertainty of 0.01 radians (about 0.57 degrees) std dev on each axis, uncorrelated.
    std_dev = 0.01
    Sigma = np.diag([std_dev**2, std_dev**2, std_dev**2])

    print("The post-reset covariance matrix Σ' is computed from the pre-reset covariance Σ using the formula:")
    print("Σ' = G * Σ * G^T\n")

    print(f"The attitude deviation vector δ to be reset is:\nδ = {delta}\n")

    print(f"The pre-reset covariance matrix Σ is:\nΣ =\n{Sigma}\n")

    # --- Calculation ---
    # The transformation of the covariance matrix is given by the formula:
    # Sigma_prime = G * Sigma * G_transpose
    # where G = exp(hat(delta)) * inv(J_l(delta))

    # 1. Compute the rotation matrix from delta
    R = exp_map(delta)

    # 2. Compute the inverse of the left Jacobian
    inv_Jl = inv_left_jacobian(delta)

    # 3. Compute the transformation matrix G
    G = R @ inv_Jl
    G_T = G.T

    # 4. Compute the post-reset covariance matrix Sigma_prime
    Sigma_prime = G @ Sigma @ G_T

    # --- Output Results ---
    print("The transformation matrix G is computed as G = exp(hat(δ)) * Jₗ(δ)⁻¹.")
    print(f"The rotation matrix exp(hat(δ)) is:\nR =\n{R}\n")
    print(f"The inverse left Jacobian Jₗ(δ)⁻¹ is:\nJₗ(δ)⁻¹ =\n{inv_Jl}\n")

    # Final equation with numbers
    print("The final equation is Σ' = G * Σ * G^T:\n")
    print("G =\n", G, "\n")
    print("Σ =\n", Sigma, "\n")
    print("G^T =\n", G_T, "\n")
    print("Resulting in the final post-reset covariance matrix Σ':\nΣ' =\n", Sigma_prime)
    
if __name__ == '__main__':
    main()
<<<The final post-reset covariance matrix Σ' is computed as Σ' = G * Σ * G^T, where G = exp(hat(δ)) * Jₗ(δ)⁻¹. The provided Python code calculates this exact expression.>>>