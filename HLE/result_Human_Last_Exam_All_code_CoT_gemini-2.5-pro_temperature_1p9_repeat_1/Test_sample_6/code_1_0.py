import numpy as np
from scipy.linalg import eigh

def solve_and_count_eigenvalues():
    """
    Solves for the Kaluza-Klein masses and counts how many are below a threshold.
    """
    # 1. Setup constants and the grid for discretization
    N = 2000  # Number of grid points for accuracy
    L = 2 * np.pi
    h = L / N
    x = np.linspace(0, L, N, endpoint=False)

    # 2. Define the warp factor and the functions for the potential V(x) and weight W(x)
    def A(t):
        return np.sin(t) + 4 * np.cos(t)

    def V(t):
        A_prime = np.cos(t) - 4 * np.sin(t)
        A_double_prime = -np.sin(t) - 4 * np.cos(t)  # Note: A'' = -A
        return 1.5 * A_double_prime + (9.0 / 4.0) * A_prime**2

    def W(t):
        return np.exp(2 * A(t))

    # Evaluate V and W on the grid points
    V_vals = V(x)
    W_vals = W(x)

    # 3. Construct the matrices H and W for the generalized eigenvalue problem H*chi = m^2*W*chi
    
    # W_mat is a diagonal matrix with the values of the weight function
    W_mat = np.diag(W_vals)

    # H_mat is the matrix representing the operator (-d^2/dx^2 + V(x)) with periodic boundary conditions
    H_mat = np.zeros((N, N))

    # Main diagonal of H_mat
    diag_vals = 2 / h**2 + V_vals
    np.fill_diagonal(H_mat, diag_vals)

    # Off-diagonals of H_mat
    off_diag_val = -1 / h**2
    np.fill_diagonal(H_mat[1:, :-1], off_diag_val)
    np.fill_diagonal(H_mat[:-1, 1:], off_diag_val)
    
    # Corner elements to enforce periodic boundary conditions
    H_mat[0, N - 1] = off_diag_val
    H_mat[N - 1, 0] = off_diag_val

    # 4. Solve the generalized eigenvalue problem. eigh returns sorted eigenvalues.
    eigenvalues, _ = eigh(H_mat, W_mat)
    
    m_squared_values = eigenvalues

    # 5. Count how many eigenvalues are below the given threshold
    threshold = 14
    count = np.sum(m_squared_values < threshold)
    
    print(f"The number of spin 2 KK mode eigenvalues (m^2) below {threshold} is:")
    print(count)

solve_and_count_eigenvalues()
<<<16>>>