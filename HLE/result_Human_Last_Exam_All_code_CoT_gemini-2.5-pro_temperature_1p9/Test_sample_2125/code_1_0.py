import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import interp1d
from scipy.optimize import brentq, root_scalar

def get_observables(alpha, x_grid):
    """
    Solves the time-independent Schrödinger equation for a given alpha
    and returns E2 and the value of psi2 at x=alpha.
    """
    if alpha <= 0:
        return np.nan, np.nan

    dx = x_grid[1] - x_grid[0]
    N = len(x_grid)

    # The potential V(x)
    V = (0.5 * alpha**2 - 3.5) * x_grid**2 - alpha * x_grid**4 + 0.5 * x_grid**6

    # Construct the Hamiltonian matrix H = T + V
    # Kinetic part T uses a three-point finite difference stencil for the second derivative.
    # H = [-1/(2*dx^2) * (psi_{i+1} - 2*psi_i + psi_{i-1})] + V_i * psi_i
    H_diag = 1.0 / dx**2 + V
    H_offdiag = -0.5 / dx**2 * np.ones(N - 1)
    
    # Using eigh_tridiagonal for performance on a tridiagonal matrix
    # We need the 3 lowest eigenvalues (n=0, 1, 2)
    try:
        eigenvalues, eigenvectors = eigh_tridiagonal(H_diag, H_offdiag, select='i', select_range=(0, 2))
    except ImportError: # For older scipy versions that may not have eigh_tridiagonal
        H = np.diag(H_diag) + np.diag(H_offdiag, k=1) + np.diag(H_offdiag, k=-1)
        eigenvalues, eigenvectors = eigh(H, subset_by_index=[0, 2])
        
    # The second excited state is the 3rd state (index 2)
    E2 = eigenvalues[2]
    psi2_on_grid = eigenvectors[:, 2]

    # Normalize the wavefunction
    psi2_on_grid /= np.sqrt(np.sum(psi2_on_grid**2) * dx)
    
    # Ensure consistent phase (e.g., positive at the center for an even function)
    center_index = N // 2
    if psi2_on_grid[center_index] < 0:
        psi2_on_grid *= -1

    # Interpolate to find psi2 at x=alpha
    if alpha >= x_grid[-1] or alpha <= x_grid[0]:
        psi2_at_alpha = 0.0
    else:
        interp_func = interp1d(x_grid, psi2_on_grid, kind='cubic')
        psi2_at_alpha = interp_func(alpha).item()

    return E2, psi2_at_alpha

def find_all_roots(func, scan_range, n_scan=200):
    """
    Finds all roots of a function within a given range.
    """
    roots = []
    alpha_grid = np.linspace(scan_range[0], scan_range[1], n_scan)
    func_vals = np.array([func(a) for a in alpha_grid])

    for i in range(len(alpha_grid) - 1):
        if np.sign(func_vals[i]) != np.sign(func_vals[i+1]):
            a1, a2 = alpha_grid[i], alpha_grid[i+1]
            try:
                sol = root_scalar(func, bracket=[a1, a2], method='brentq')
                if sol.converged:
                    roots.append(sol.root)
            except (ValueError, RuntimeError):
                continue
    return roots

def main():
    """
    Main function to perform the calculation and print the result.
    """
    # Grid setup for the numerical solver
    x_min = -10.0
    x_max = 10.0
    N_grid = 4001 # Use an odd number to have a grid point at x=0
    x_grid = np.linspace(x_min, x_max, N_grid)

    # Define the functions whose roots we want to find
    def func_E2(alpha):
        E2, _ = get_observables(alpha, x_grid)
        return E2

    def func_psi2(alpha):
        _, psi2_at_alpha = get_observables(alpha, x_grid)
        return psi2_at_alpha
        
    # Define search range for alpha
    alpha_scan_range = (0.1, 6.0)

    # Find roots for E2(alpha) = 0
    roots_E2 = find_all_roots(func_E2, alpha_scan_range)
    
    # Find roots for psi2(alpha; alpha) = 0
    roots_psi2 = find_all_roots(func_psi2, alpha_scan_range)

    # Combine all found roots
    all_roots = roots_E2 + roots_psi2
    
    if not all_roots:
        print("No root found for F(alpha) = 0 in the specified range.")
    else:
        alpha_0 = max(all_roots)
        # Final answer statement as requested by the problem description format
        print(f"The largest value alpha_0 such that F(alpha_0) = 0 is alpha_0 = {alpha_0:.10f}")
        print(f"<<<{alpha_0:.10f}>>>")

if __name__ == '__main__':
    main()
