import numpy as np

def skew_symmetric(v):
    """
    Creates a skew-symmetric matrix from a 3-element vector.
    hat(v)
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def left_jacobian_so3(delta):
    """
    Computes the left Jacobian of SO(3) for a given rotation vector delta.
    """
    theta = np.linalg.norm(delta)
    delta_skew = skew_symmetric(delta)
    
    # Use Taylor series expansion for small theta to avoid division by zero
    if np.isclose(theta, 0.0):
        # J_l(0) = I + 1/2 * hat(delta) + 1/6 * hat(delta)^2
        A = 1.0 / 2.0
        B = 1.0 / 6.0
    else:
        theta2 = theta * theta
        theta3 = theta2 * theta
        A = (1 - np.cos(theta)) / theta2
        B = (theta - np.sin(theta)) / theta3
        
    I = np.identity(3)
    J_l = I + A * delta_skew + B * (delta_skew @ delta_skew)
    return J_l

def main():
    """
    Main function to demonstrate the covariance reset calculation.
    """
    # --- Inputs ---
    # Pre-reset attitude deviation vector (the filter's estimate)
    delta = np.array([0.1, -0.1, 0.2])
    
    # Pre-reset covariance matrix of the attitude deviation
    # (assuming variances of (0.1 rad)^2 on x,y and (0.15 rad)^2 on z)
    Sigma = np.diag([0.1**2, 0.1**2, 0.15**2])
    
    # --- Calculation ---
    # 1. Compute the left Jacobian matrix
    J_l = left_jacobian_so3(delta)
    
    # 2. Compute the post-reset covariance using the formula
    Sigma_prime = J_l @ Sigma @ J_l.T
    
    # --- Output ---
    print("This script calculates the post-reset covariance Σ' for an attitude deviation.")
    print("The transformation is: Σ' = J_l(δ) * Σ * J_l(δ)^T\n")
    
    np.set_printoptions(precision=6, suppress=True)
    
    print("Given the pre-reset deviation vector δ:")
    print(delta)
    print("\nAnd the pre-reset covariance Σ:")
    print(Sigma)
    
    print("\n1. The computed Left Jacobian matrix J_l(δ) is:")
    print(J_l)
    
    print("\n2. The post-reset covariance Σ' is calculated as J_l * Σ * J_l^T:")
    print(Sigma_prime)
    
    print("\n--- Equation with all numbers ---")
    print("\nΣ' = \n", J_l, "\n\n  *  \n\n", Sigma, "\n\n  *  \n\n", J_l.T)
    print("\n\nResulting Σ' = \n", Sigma_prime)


if __name__ == "__main__":
    main()