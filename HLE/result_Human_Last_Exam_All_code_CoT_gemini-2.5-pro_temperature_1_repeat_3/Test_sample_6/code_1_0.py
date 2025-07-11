import numpy as np
from scipy.linalg import eigh

def solve_kk_masses():
    """
    This function calculates the Kaluza-Klein mass spectrum for a given
    5D gravity model and counts the number of modes below a threshold.
    """
    # 1. Set up the numerical parameters and the grid
    N = 500  # Number of grid points for discretization
    L = 2 * np.pi  # The length of the interval for the coordinate x
    dx = L / N
    x = np.linspace(0, L, N, endpoint=False)

    # 2. Define the functions from the physical setup
    # The warp factor function A(x)
    A = np.sin(x) + 4 * np.cos(x)
    # The functions p(x) and w(x) for the Sturm-Liouville problem
    p_func = np.exp(A)
    w_func = np.exp(3 * A)

    # 3. Construct the matrices for the generalized eigenvalue problem M*v = lambda*W*v
    
    # W is a diagonal matrix with the values of w_func on the diagonal
    W = np.diag(w_func)

    # M is a symmetric matrix representing the differential operator.
    # We evaluate p(x) at the midpoints of the grid for better accuracy.
    x_half_step = x + dx / 2
    A_half_step = np.sin(x_half_step) + 4 * np.cos(x_half_step)
    p_half_step = np.exp(A_half_step)

    # Main diagonal of M
    M_diag = (np.roll(p_half_step, 1) + p_half_step) / dx**2
    # Off-diagonal elements of M
    M_offdiag = -p_half_step / dx**2
    
    # Assemble the matrix M
    M = np.diag(M_diag) + np.diag(M_offdiag[:-1], k=1) + np.diag(M_offdiag[:-1], k=-1)
    
    # Add corner elements for periodic boundary conditions
    M[0, N - 1] = M_offdiag[-1]
    M[N - 1, 0] = M_offdiag[-1]

    # 4. Solve the eigenvalue problem
    # eigh is used for symmetric/Hermitian matrices and is numerically stable.
    eigenvalues, eigenvectors = eigh(M, W)
    
    # Sort eigenvalues for clarity
    eigenvalues.sort()

    # 5. Count the number of eigenvalues below the threshold
    threshold = 14.0
    modes_below_threshold = eigenvalues[eigenvalues < threshold]
    count = len(modes_below_threshold)

    # 6. Print the results
    print("The eigenvalues (masses-squared) below 14 are:")
    # We print each eigenvalue as requested by the prompt
    print(modes_below_threshold)
    print(f"\nThe total number of eigenvalues below {threshold} is: {count}")
    
    return count

# Run the calculation and store the final answer
final_count = solve_kk_masses()

# The final answer in the required format
# print(f"\n<<<{final_count}>>>") # This is for internal formatting, the user just sees the result.

if __name__ == '__main__':
    # This block is for executing the code. The final answer is wrapped in <<<>>>
    # The printed output from the function is already sufficient for the user.
    # The line below is a placeholder for the final answer extraction logic.
    pass