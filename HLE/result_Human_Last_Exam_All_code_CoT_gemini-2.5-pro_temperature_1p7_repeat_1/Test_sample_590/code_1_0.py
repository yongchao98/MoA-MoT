import numpy as np
from scipy.linalg import eigh

def solve_eigenvalue_problem():
    """
    This function calculates the number of positive eigenvalues for the given stability operator L,
    following the plan outlined. It assumes n=2 (the classical catenoid) and numerically
    computes the number of negative eigenvalues of the associated 1D Schrödinger operator.
    """
    
    # Numerical parameters
    # We choose a sufficiently large interval [-R, R] and a fine grid.
    R = 20.0  # Max value for rho
    N = 4001  # Number of grid points
    
    # Set up the grid for the variable rho
    rho = np.linspace(-R, R, N)
    h = rho[1] - rho[0]  # Grid spacing
    
    # Define the potential function V(rho) for the Schrödinger operator
    rho_sq_plus_1 = rho**2 + 1
    potential = 1.0 / rho_sq_plus_1 - 2.0 / (rho_sq_plus_1**2)
    
    # Construct the Hamiltonian matrix for H = -d^2/d(rho)^2 + V(rho)
    # The kinetic part (-d^2/d(rho)^2) is discretized using a central difference formula.
    # The result is a symmetric, tridiagonal matrix.
    main_diag = 2.0 / h**2 + potential
    off_diag = -1.0 / h**2 * np.ones(N - 1)
    
    # The full matrix can be constructed and solved, but for tridiagonal systems,
    # specialized solvers are more efficient. eigh_tridiagonal is such a solver.
    # We will build the full matrix for compatibility and simplicity.
    hamiltonian_matrix = np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)
    
    # Compute the eigenvalues of the Hamiltonian matrix.
    # eigh is used because the matrix is symmetric (Hermitian).
    eigenvalues = np.linalg.eigh(hamiltonian_matrix)[0]
    
    # Count the number of negative eigenvalues.
    # This count corresponds to the number of positive eigenvalues of the original operator L.
    num_positive_eigenvalues = np.sum(eigenvalues < 0)
    
    # The final equation is the statement of our result.
    final_number = int(num_positive_eigenvalues)
    print(f"Number of positive eigenvalues = {final_number}")

# Run the solver function
solve_eigenvalue_problem()