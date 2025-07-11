import numpy as np

def count_kk_modes():
    """
    This script calculates the number of spin-2 Kaluza-Klein mode masses-squared
    (eigenvalues) below a certain threshold in a warped 5D compactification.

    The problem is mapped to solving a 1D Schrödinger equation numerically.
    """
    
    # 1. Define the warp factor and the resulting potential
    
    # The warp factor A(x) and its derivatives
    def A_func(x):
        return np.sin(x) + 4 * np.cos(x)

    def A_prime(x):
        return np.cos(x) - 4 * np.sin(x)

    def A_double_prime(x):
        return -np.sin(x) - 4 * np.cos(x)

    # The potential V(x) in the Schrödinger equation
    def potential(x):
        return (9.0 / 4.0) * A_prime(x)**2 + (3.0 / 2.0) * A_double_prime(x)

    # 2. Set up the numerical grid and parameters
    
    # Discretization parameters for the finite difference method
    N = 2000  # Number of grid points (higher for more accuracy)
    L = 2 * np.pi  # The domain is a circle of circumference 2*pi
    h = L / N  # Grid spacing
    
    # The grid of x values, excluding the endpoint for periodicity
    x_grid = np.linspace(0, L, N, endpoint=False)
    
    # 3. Construct the Hamiltonian matrix for the operator H = -d^2/dx^2 + V(x)
    
    # Kinetic part: The matrix for -d^2/dx^2 with periodic boundary conditions
    K = np.zeros((N, N))
    diag_val = 2.0 / h**2
    off_diag_val = -1.0 / h**2
    
    np.fill_diagonal(K, diag_val)
    np.fill_diagonal(K[1:], off_diag_val)
    np.fill_diagonal(K[:, 1:], off_diag_val)
    
    # Apply periodic boundary conditions at the corners
    K[0, N - 1] = off_diag_val
    K[N - 1, 0] = off_diag_val

    # Potential part: A diagonal matrix with the potential values on the diagonal
    V_matrix = np.diag(potential(x_grid))

    # The full Hamiltonian matrix
    H = K + V_matrix

    # 4. Find eigenvalues and count those below the threshold
    
    # Calculate the eigenvalues of the Hermitian matrix H
    eigenvalues = np.linalg.eigvalsh(H)
    
    # Define the threshold for the squared mass
    threshold = 14.0
    
    # Find all eigenvalues below the threshold
    eigenvalues_below_threshold = eigenvalues[eigenvalues < threshold]
    
    # Count them (degenerate eigenvalues are counted by their multiplicity)
    count = len(eigenvalues_below_threshold)
    
    # 5. Output the results
    
    print("The eigenvalues (squared masses m^2) found below the threshold of 14 are:")
    # We sort them for a clear presentation
    for E in sorted(eigenvalues_below_threshold):
        # The result includes one eigenvalue for m^2=0, and two pairs of degenerate eigenvalues.
        print(f"{E:.4f}")

    print(f"\nThe total count of eigenvalues below {threshold} is:")
    print(count)


# Execute the function to solve the task
count_kk_modes()