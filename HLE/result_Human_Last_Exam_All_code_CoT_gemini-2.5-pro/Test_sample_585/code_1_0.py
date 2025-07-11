import numpy as np

def skew_symmetric(v):
    """
    Computes the skew-symmetric matrix for a 3-element vector v.
    hat(v) = [ 0 -v3  v2]
             [ v3  0 -v1]
             [-v2  v1  0]
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def right_jacobian_so3(delta):
    """
    Computes the right Jacobian of SO(3) for a given rotation vector delta.
    """
    norm_delta = np.linalg.norm(delta)
    
    # Handle the case of a very small rotation to avoid division by zero
    if norm_delta < 1e-9:
        return np.identity(3)
        
    I = np.identity(3)
    hat_delta = skew_symmetric(delta)
    hat_delta_sq = hat_delta @ hat_delta
    
    norm_delta_sq = norm_delta**2
    norm_delta_cub = norm_delta**3
    
    cos_norm = np.cos(norm_delta)
    sin_norm = np.sin(norm_delta)
    
    # Coefficients for the Jacobian formula
    A = (1 - cos_norm) / norm_delta_sq
    B = (norm_delta - sin_norm) / norm_delta_cub
    
    # The right Jacobian of SO(3)
    Jr = I - A * hat_delta + B * hat_delta_sq
    return Jr

# --- Main execution ---

# Define a sample pre-reset attitude deviation vector (the filter's estimate)
# This represents a rotation of about 0.24 rad (14 degrees)
delta = np.array([0.1, -0.1, 0.2])

# Define a sample pre-reset covariance matrix for the deviation vector
# This assumes small uncertainty, with variances of 0.01^2 and no correlation
Sigma = np.diag([0.01**2, 0.01**2, 0.01**2])

# 1. Compute the right Jacobian J_r(delta)
Jr = right_jacobian_so3(delta)

# 2. Compute the post-reset covariance Sigma'
Sigma_prime = Jr @ Sigma @ Jr.T

# Print all the inputs and results
print("--- Input Values ---")
print("Pre-reset deviation (delta):\n", delta)
print("\nPre-reset covariance (Sigma):\n", Sigma)

print("\n\n--- Calculation Results ---")
print("Right Jacobian of SO(3) (Jr):\n", Jr)
print("\nPost-reset covariance (Sigma'):\n", Sigma_prime)
