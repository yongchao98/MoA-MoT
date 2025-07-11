import numpy as np
import scipy.linalg as linalg
from scipy.special import binom

def solve_eigenvalues():
    """
    Solves for the number of positive eigenvalues for the given stability operator L.
    The primary case for a catenoid corresponds to n=2, which we will use.
    """
    n = 2  # Dimension of the catenoid, n=2 for a surface in R^3.

    # Numerical parameters for the solver
    R = 20.0  # Domain for rho is [-R, R]
    N = 2000  # Number of grid points

    rho = np.linspace(-R, R, N)
    h = rho[1] - rho[0]  # Grid spacing

    # The potential is positive only for k < n-1.
    # For n=2, this is only k=0.
    # The multiplicity of k=0 spherical harmonic is 1.
    
    total_positive_eigenvalues = 0
    k = 0 # We only need to check the k=0 mode for n=2

    # Define helper functions based on the problem description
    # We assume |F_rho| has a typo and should contain abs(rho).
    # Small epsilon to avoid division by zero at rho=0
    eps = 1e-9
    
    def japanese_bracket(rho):
        return np.sqrt(rho**2 + 1)

    def F_rho_abs(rho_val, n_val):
        rho_val = np.abs(rho_val)
        jb = japanese_bracket(rho_val)
        jb_pow = jb**(2 * (n_val - 1))
        # Handle the case where jb_pow is close to 1
        denom = np.sqrt(np.maximum(jb_pow - 1, eps))
        return (rho_val * jb**(n_val - 2)) / denom

    # Coefficients for the Sturm-Liouville problem (p(rho)u')' + V(rho)w(rho)u = lambda w(rho)u
    jb = japanese_bracket(rho)
    
    # We need p(rho) = <rho>^(n-1) / |F_rho| and w(rho) = <rho>^(n-1) * |F_rho|
    frho = F_rho_abs(rho, n)
    p_rho = jb**(n-1) / frho
    w_rho = jb**(n-1) * frho
    
    # Potential V_k(rho)
    V_k = -k*(k+n-2)/jb**2 + n*(n-1)/jb**(2*n)

    # Assemble the matrices for the generalized eigenvalue problem Ax = lambda Bx
    # using a finite element method with piecewise linear basis functions (hat functions).
    # This leads to symmetric tridiagonal matrices A and B.

    # Matrix for the kinetic term: integral(-p(rho) * phi_i' * phi_j')
    A_kin = np.zeros((N - 2, N - 2))
    p_mid = 0.5 * (p_rho[:-1] + p_rho[1:]) # p at midpoints
    diag_kin = -(p_mid[:-1] + p_mid[1:]) / h
    offdiag_kin = p_mid[1:-1] / h
    np.fill_diagonal(A_kin, diag_kin)
    np.fill_diagonal(A_kin[1:], offdiag_kin)
    np.fill_diagonal(A_kin[:, 1:], offdiag_kin)

    # Matrix for the potential term: integral(V_k(rho) * w(rho) * phi_i * phi_j)
    A_pot = np.zeros((N-2, N-2))
    Vw = V_k * w_rho
    Vw_int_diag = (2.0/3.0) * h * Vw[1:-1]
    Vw_int_offdiag = (1.0/6.0) * h * (Vw[1:-2] + Vw[2:-1])
    np.fill_diagonal(A_pot, Vw_int_diag)
    np.fill_diagonal(A_pot[1:], Vw_int_offdiag)
    np.fill_diagonal(A_pot[:, 1:], Vw_int_offdiag)
    
    # Total A matrix
    A = A_kin + A_pot

    # Mass matrix B: integral(w(rho) * phi_i * phi_j)
    B = np.zeros((N - 2, N - 2))
    w_int_diag = (2.0/3.0) * h * w_rho[1:-1]
    w_int_offdiag = (1.0/6.0) * h * (w_rho[1:-2] + w_rho[2:-1])
    np.fill_diagonal(B, w_int_diag)
    np.fill_diagonal(B[1:], w_int_offdiag)
    np.fill_diagonal(B[:, 1:], w_int_offdiag)

    # Solve the generalized eigenvalue problem
    # We only need the eigenvalues, so we can use eigvalsh
    eigenvalues = linalg.eigvalsh(A, B)

    # Count positive eigenvalues
    num_positive = np.sum(eigenvalues > 0)

    # For n=2, multiplicity of k=0 is 1
    total_positive_eigenvalues += num_positive
            
    print(f"For n={n}, the number of positive eigenvalues is determined numerically.")
    print(f"The calculation for k=0 gives {num_positive} positive eigenvalue(s).")
    print(f"Modes with k > 0 do not contribute positive eigenvalues for n=2.")
    print(f"Total number of positive eigenvalues is: {total_positive_eigenvalues}")

solve_eigenvalues()
