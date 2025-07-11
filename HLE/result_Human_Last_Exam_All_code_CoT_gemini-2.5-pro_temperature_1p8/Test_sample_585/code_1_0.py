import numpy as np

def skew_symmetric(v):
    """
    Creates a skew-symmetric matrix from a 3-element vector.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def left_jacobian_so3(phi):
    """
    Computes the left Jacobian of SO(3).

    J_l(phi) = I + (1 - cos(theta))/theta^2 * hat(phi) + (theta - sin(theta))/theta^3 * hat(phi)^2

    Args:
        phi: A 3-element numpy array representing the rotation vector.

    Returns:
        A 3x3 numpy array for the left Jacobian.
    """
    theta = np.linalg.norm(phi)
    phi_skew = skew_symmetric(phi)
    
    # Use small angle approximation to avoid division by zero
    if np.isclose(theta, 0.0):
        # Taylor series expansion for the coefficients
        A = 0.5 - theta**2 / 24.0
        B = 1.0/6.0 - theta**2 / 120.0
    else:
        theta2 = theta * theta
        theta3 = theta2 * theta
        A = (1 - np.cos(theta)) / theta2
        B = (theta - np.sin(theta)) / theta3
        
    I = np.identity(3)
    J_l = I + A * phi_skew + B * np.dot(phi_skew, phi_skew)
    
    return J_l

def main():
    """
    Demonstrates the covariance update after an attitude reset.
    """
    # Let's assume the attitude deviation `delta` to be reset is:
    # (this is the mean of the error distribution from the Kalman filter)
    delta = np.array([0.1, -0.05, 0.02])

    # Let's assume the pre-reset covariance `Sigma` is:
    # (a diagonal matrix for simplicity, representing uncertainty in each axis)
    Sigma = np.diag([0.01**2, 0.01**2, 0.02**2])

    print("--- Attitude Reset Covariance Propagation ---")
    print("\nPre-reset state:")
    print(f"Attitude deviation (delta) = {delta}")
    print("\nPre-reset covariance (Sigma):")
    print(Sigma)
    
    # 1. Compute the Left Jacobian of SO(3) for `delta`
    J_l = left_jacobian_so3(delta)

    # 2. Compute the post-reset covariance using the formula: Sigma' = J_l * Sigma * J_l^T
    Sigma_prime = J_l @ Sigma @ J_l.T

    print("\nTransformation Jacobian (J_l(delta)):")
    print(J_l)

    print("\nPost-reset covariance (Sigma'):")
    print(Sigma_prime)

    print("\nThe final equation is Sigma' = J_l(delta) * Sigma * J_l(delta)^T")

if __name__ == '__main__':
    main()
>>>