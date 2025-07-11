import numpy as np
from scipy.linalg import eigh

def solve_kk_masses():
    """
    This function calculates the eigenvalues for the spin-2 Kaluza-Klein modes
    and counts how many are below a specified threshold.
    """
    # Define the warp factor and derived functions p(x) and w(x)
    def A(x):
        return np.sin(x) + 4 * np.cos(x)

    def p(x):
        return np.exp(3 * A(x))

    def w(x):
        return np.exp(A(x))

    # Set up the numerical discretization
    N = 2000  # Number of grid points for high accuracy
    h = 2 * np.pi / N  # Grid spacing

    # Create the grid for x in [0, 2*pi)
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)

    # Construct the matrices L and W for the generalized eigenvalue problem L*v = m^2*W*v
    
    # W is a diagonal matrix with the weight function values
    w_vals = w(x)
    W = np.diag(w_vals)

    # L is a symmetric matrix representing the differential operator -(p*psi')'
    # We use a finite difference scheme on a staggered grid for p(x)
    x_half = x + h / 2
    p_half_vals = p(x_half)

    # The main diagonal of L
    diag_L = (p_half_vals + np.roll(p_half_vals, 1)) / h**2
    # The off-diagonals of L
    off_diag_L = -p_half_vals / h**2

    # Assemble the matrix L with periodic boundary conditions
    L = np.diag(diag_L) + np.diag(off_diag_L[:-1], k=1) + np.diag(off_diag_L[:-1], k=-1)
    L[0, N - 1] = off_diag_L[-1]
    L[N - 1, 0] = off_diag_L[-1]

    # Solve the generalized eigenvalue problem for the eigenvalues (m^2)
    # eigh is used for real symmetric/Hermitian matrices and returns sorted eigenvalues.
    eigenvalues, eigenvectors = eigh(L, W)

    # Define the threshold and find the eigenvalues below it
    threshold = 14
    below_threshold_e_vals = eigenvalues[eigenvalues < threshold]
    count = len(below_threshold_e_vals)

    # Output the results as requested
    print(f"The eigenvalues (m^2) below the threshold of {threshold} are:")
    # The instruction "output each number in the final equation" is interpreted
    # as printing the eigenvalues that contribute to the final count.
    for val in below_threshold_e_vals:
        print(f"{val:.6f}")

    print(f"\nThe total number of eigenvalues below {threshold} is:")
    print(count)
    
    return count

# Run the calculation and store the final answer
final_count = solve_kk_masses()
# The final answer is provided in the specified format below.
# <<<37>>>