import numpy as np

def solve_kk_mass_spectrum():
    """
    This function calculates the number of Kaluza-Klein eigenvalues below a
    certain threshold for a given 5D gravitational theory.
    """
    # Step 1: Define the problem parameters and the grid.
    N = 2000  # Number of grid points for discretization
    L = 2 * np.pi  # The length of the interval for x
    h = L / N  # Grid spacing
    x = np.linspace(0, L, N, endpoint=False)  # Grid points from 0 to 2*pi

    # Step 2: Define the warp factor A(x) and its derivatives.
    # The warp factor is given by A(x) = sin(x) + 4*cos(x).
    A = np.sin(x) + 4 * np.cos(x)
    # The first derivative is A'(x) = cos(x) - 4*sin(x).
    A_prime = np.cos(x) - 4 * np.sin(x)
    # The second derivative is A''(x) = -sin(x) - 4*cos(x).
    A_double_prime = -np.sin(x) - 4 * np.cos(x)

    # Step 3: Construct the potential V(x).
    # The potential for spin-2 modes is V(x) = (3/2)A''(x) + (9/4)(A'(x))^2.
    V = (3 / 2) * A_double_prime + (9 / 4) * (A_prime**2)

    # Step 4: Construct the Hamiltonian matrix for the discretized operator.
    # The operator is H = -d^2/dx^2 + V(x).
    # We use a centered finite difference scheme for the second derivative
    # with periodic boundary conditions.
    
    # Diagonal part of the Hamiltonian
    H_diag = V + 2 / h**2
    H = np.diag(H_diag)

    # Off-diagonal parts
    off_diag_val = -1 / h**2
    H += np.diag(off_diag_val * np.ones(N - 1), k=1)
    H += np.diag(off_diag_val * np.ones(N - 1), k=-1)

    # Corner elements for periodic boundary conditions
    H[0, N - 1] = off_diag_val
    H[N - 1, 0] = off_diag_val

    # Step 5: Calculate the eigenvalues of the Hamiltonian matrix.
    # We use eigvalsh since the matrix is symmetric (Hermitian).
    eigenvalues = np.linalg.eigvalsh(H)

    # Step 6: Count the number of eigenvalues below the threshold of 14.
    threshold = 14
    count = np.sum(eigenvalues < threshold)

    # Output the results
    print("The eigenvalue problem to be solved for the squared mass m^2 is:")
    print("-f''(x) + V(x)f(x) = m^2 f(x)")
    print("where the warp factor A(x) and the potential V(x) are:")
    print("A(x) = sin(x) + 4*cos(x)")
    print("V(x) = (3/2)*A''(x) + (9/4)*(A'(x))^2")
    print("\nNumerically solving for the eigenvalues (m^2)...")
    
    # Print the lowest few eigenvalues to show the structure
    sorted_eigenvalues = np.sort(eigenvalues)
    print("\nThe lowest computed eigenvalues are:")
    for i in range(min(10, len(sorted_eigenvalues))):
        print(f"m_{i}^2 = {sorted_eigenvalues[i]:.4f}")

    print(f"\nCounting the number of eigenvalues below {threshold}...")
    print(f"The number of eigenvalues below {threshold} is: {count}")

if __name__ == '__main__':
    solve_kk_mass_spectrum()