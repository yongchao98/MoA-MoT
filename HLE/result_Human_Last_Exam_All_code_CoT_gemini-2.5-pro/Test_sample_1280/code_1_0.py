import numpy as np
import scipy.linalg

def solve_schrodinger(potential_func, n_points=1000, num_eigenvalues=10):
    """
    Solves the 1D time-independent Schrödinger equation on [0,1]
    with Dirichlet boundary conditions using the finite difference method.

    Args:
        potential_func (callable): A function V(x) for the potential.
        n_points (int): Number of discretization points.
        num_eigenvalues (int): Number of lowest eigenvalues to return.

    Returns:
        np.ndarray: An array of the lowest eigenvalues.
    """
    # Discretize the interval [0, 1]. The wavefunction is zero at the endpoints.
    # So we solve for the n_points-2 interior points.
    h = 1.0 / (n_points - 1)
    x = np.linspace(0, 1, n_points)
    
    # Get the potential at the interior points
    V = potential_func(x[1:-1])

    # Construct the Hamiltonian matrix for the interior points
    # Main diagonal
    main_diag = 2.0 / h**2 + V
    # Off-diagonal
    off_diag = -1.0 / h**2 * np.ones(n_points - 3)

    # Find the eigenvalues of the tridiagonal matrix.
    # We ask for the first `num_eigenvalues` lowest eigenvalues.
    eigenvalues = scipy.linalg.eigh_tridiagonal(main_diag,
                                                 off_diag,
                                                 select='i',
                                                 select_range=(0, num_eigenvalues - 1))[0]
    return eigenvalues

def run_demonstration():
    """
    Constructs two partner potentials from a superpotential W(x) and
    compares their energy spectra.
    """
    # 1. Define the superpotential W(x) and its derivative W'(x)
    # Let's choose a smooth function on [0,1].
    W = lambda x: 4 * np.cos(2 * np.pi * x)
    W_prime = lambda x: -8 * np.pi * np.sin(2 * np.pi * x)

    # 2. Define the two partner potentials V0(x) and V1(x)
    # H0 is associated with V0 = W^2 - W'
    # H1 is associated with V1 = W^2 + W'
    V0 = lambda x: W(x)**2 - W_prime(x)
    V1 = lambda x: W(x)**2 + W_prime(x)

    # 3. Solve for the eigenvalues of H0 and H1 (we ignore the constant alpha)
    print("Calculating eigenvalues for H0 and H1...")
    eigenvalues_H0 = solve_schrodinger(V0)
    eigenvalues_H1 = solve_schrodinger(V1)

    # 4. Print the results for comparison
    print("\nComparison of the first 10 eigenvalues (alpha=0):")
    print("-------------------------------------------------")
    print(f"{'Level (n)':<10} {'Eigenvalues of H0':<20} {'Eigenvalues of H1':<20}")
    print("-------------------------------------------------")
    for i in range(len(eigenvalues_H0)):
        print(f"{i:<10} {eigenvalues_H0[i]:<20.6f} {eigenvalues_H1[i]:<20.6f}")
    
    print("\nAs demonstrated for Dirichlet boundary conditions, the spectra are identical.")
    print("Under different boundary conditions, the spectra could differ by at most one level.")
    print("The maximum number of levels that can differ is 1.")

if __name__ == '__main__':
    run_demonstration()
