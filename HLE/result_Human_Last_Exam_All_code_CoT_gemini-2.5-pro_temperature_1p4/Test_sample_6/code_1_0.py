import numpy as np
from scipy.linalg import eigh

def solve_kk_masses():
    """
    This function calculates the eigenvalues for the Kaluza-Klein modes
    and counts how many are below a given threshold.
    """
    # 1. Setup the problem parameters
    # Using a sufficiently large number of grid points for accuracy
    N = 400
    L = 2 * np.pi
    h = L / N
    x = np.linspace(0, L, N, endpoint=False)  # Grid for the circle coordinate x

    # The threshold for the squared mass m^2
    threshold = 14

    # 2. Define the warp factor and the related function p(x)
    # The warp factor A(x)
    def A(x):
        return np.sin(x) + 4 * np.cos(x)

    # The function p(x) = w(x) = exp(4A(x)) from the Sturm-Liouville equation
    def p(x):
        return np.exp(4 * A(x))

    # Evaluate p(x) on the grid
    p_values = p(x)

    # 3. Construct the matrices for the generalized eigenvalue problem A v = lambda B v
    # B is the diagonal mass matrix from the right-hand side, representing the weight function w(x)
    B = np.diag(p_values)

    # A is the matrix representing the differential operator L = -d/dx (p(x) d/dx)
    # This uses a centered finite difference scheme that preserves the self-adjointness.
    
    # First, calculate p at midpoints p_{j+1/2}
    p_half = (p_values + np.roll(p_values, -1)) / 2.0

    # The main diagonal of matrix A is (p_{j+1/2} + p_{j-1/2}) / h^2
    main_diag = (p_half + np.roll(p_half, 1)) / h**2

    # The off-diagonals are -p_{j+1/2} / h^2
    off_diag = -p_half / h**2

    # Assemble the matrix A
    A_matrix = np.diag(main_diag) + np.diag(off_diag[:-1], 1) + np.diag(off_diag[:-1], -1)

    # Add corner elements to enforce periodic boundary conditions
    A_matrix[0, N - 1] = off_diag[-1]
    A_matrix[N - 1, 0] = off_diag[-1]

    # 4. Solve the generalized eigenvalue problem
    # scipy.linalg.eigh is used for symmetric/Hermitian matrices and returns
    # eigenvalues in ascending order.
    eigenvalues, eigenvectors = eigh(A_matrix, B)

    # 5. Filter, count, and print the results
    print("The squared masses m^2 of the Kaluza-Klein modes are the eigenvalues of the equation:")
    print("-(exp(4*(sin(x)+4*cos(x))) * psi'(x))' = m^2 * exp(4*(sin(x)+4*cos(x))) * psi(x)")
    print(f"\nWe need to count the number of eigenvalues below {threshold}.")
    print("The eigenvalues found to be below the threshold are:")

    # Filter eigenvalues below the threshold
    below_threshold_eigenvalues = eigenvalues[eigenvalues < threshold]
    
    for val in below_threshold_eigenvalues:
        print(f"{val:.4f}")
    
    count = len(below_threshold_eigenvalues)
    
    print(f"\nThe total number of eigenvalues below {threshold} is {count}.")
    
    # Return the final count in the specified format
    print(f"\n<<<{count}>>>")

if __name__ == '__main__':
    solve_kk_masses()