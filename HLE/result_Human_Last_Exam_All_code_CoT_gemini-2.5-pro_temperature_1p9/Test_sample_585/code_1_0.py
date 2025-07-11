import numpy as np

def skew(v):
    """
    Computes the skew-symmetric matrix for a 3-element vector v.
    S(v) = [[  0, -v3,  v2],
            [ v3,   0, -v1],
            [-v2,  v1,   0]]
    """
    v = v.flatten()
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def left_jacobian_so3(phi):
    """
    Computes the Left Jacobian of SO(3) for a given rotation vector phi.
    J_l(phi) = I + (1-cos(||phi||))/||phi||^2 * S(phi) + (||phi||-sin(||phi||))/||phi||^3 * S(phi)^2
    """
    phi = phi.flatten()
    theta = np.linalg.norm(phi)
    
    # Handle the case of a very small rotation to avoid division by zero
    if theta < 1e-9:
        return np.identity(3)
        
    S = skew(phi)
    S2 = S @ S
    
    theta2 = theta * theta
    theta3 = theta2 * theta
    
    # Coefficients for the Jacobian formula
    A = (1 - np.cos(theta)) / theta2
    B = (theta - np.sin(theta)) / theta3
    
    # The Left Jacobian matrix
    J_l = np.identity(3) + A * S + B * S2
    return J_l

def format_matrix(m, name):
    """Helper function to format a matrix for printing."""
    header = f"{name} ="
    # Use a generic format for floating point numbers
    matrix_str = np.array2string(m, formatter={'float_kind':lambda x: "%.6f" % x})
    return f"{header}\n{matrix_str}"
    
# --- Example Case ---
# The reset vector 'delta', which is the mean of the attitude deviation state.
# For example, a rotation of approx 0.2 rads (~11.5 deg).
delta_hat = np.array([0.1, -0.15, 0.08])

# The pre-reset covariance matrix 'Sigma'.
# Let's assume some small, non-diagonal covariance.
Sigma = np.array([
    [0.0004, 0.00005, -0.0001],
    [0.00005, 0.0005, 0.00002],
    [-0.0001, 0.00002, 0.0003]
])

# 1. Compute the Left Jacobian matrix for the reset vector
G = left_jacobian_so3(delta_hat)

# 2. Compute the post-reset covariance matrix 'Sigma_prime'
Sigma_prime = G @ Sigma @ G.T

# --- Output the results ---
print("This script demonstrates the calculation of the post-reset attitude covariance.")
print("-" * 70)
print(f"Given a reset vector delta_hat:\n{delta_hat}\n")
print(format_matrix(Sigma, "And a pre-reset covariance matrix Sigma"))
print("-" * 70)
print("The post-reset covariance Sigma_prime is computed as:\n")
print("Sigma_prime = G * Sigma * G^T\n")
print("Where G is the Left Jacobian of SO(3) evaluated at delta_hat.\n")

print(format_matrix(G, "Computed Left Jacobian, G"))
print("\n*\n")
print(format_matrix(Sigma, "Sigma (pre-reset)"))
print("\n*\n")
print(format_matrix(G.T, "G Transposed"))
print("\n=\n")
print(format_matrix(Sigma_prime, "Sigma_prime (post-reset)"))
print("-" * 70)
