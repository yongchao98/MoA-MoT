import numpy as np
import scipy.special

def compute_legendre_roots_weights(n):
    # Use scipy to compute the roots and weights for the Legendre polynomial
    roots, weights = np.polynomial.legendre.leggauss(n)
    
    # Filter out only the non-negative roots and corresponding weights
    non_negative_indices = roots >= 0
    non_negative_roots = roots[non_negative_indices]
    non_negative_weights = weights[non_negative_indices]
    
    return non_negative_roots.tolist(), non_negative_weights.tolist()

def main_solution(n):
    roots, weights = compute_legendre_roots_weights(n)
    return {"roots": roots, "weights": weights}

# Test the function with n = 9
result = main_solution(9)
print(result)