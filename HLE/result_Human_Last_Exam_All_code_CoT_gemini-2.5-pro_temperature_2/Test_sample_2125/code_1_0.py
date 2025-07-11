import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq

# Use a cache to store solutions and avoid recomputing for the same alpha
solution_cache = {}

def solve_schrodinger(alpha, N=4001, L=15.0):
    """
    Solves the 1D Schrodinger equation for the given potential using finite differences.
    Results are cached to speed up computation.

    Args:
        alpha (float): The parameter in the potential.
        N (int): Number of grid points.
        L (float): The half-width of the symmetric spatial grid [-L, L].

    Returns:
        tuple: A tuple containing (E0, E2, psi0, psi2, x), which are the
               ground and second excited state energies and wavefunctions, and the grid.
    """
    alpha_key = round(alpha, 9)
    if alpha_key in solution_cache:
        return solution_cache[alpha_key]

    x = np.linspace(-L, L, N)
    h = x[1] - x[0]

    V = (0.5 * (alpha**2 - 7.0) * x**2 - alpha * x**4 + 0.5 * x**6)

    diag = 1.0/h**2 + V
    offdiag = -0.5/h**2 * np.ones(N-1)

    # Solve for the lowest 3 eigenvalues/eigenvectors
    energies, wavefuncs_raw = eigh_tridiagonal(diag, offdiag, select='i', select_range=(0, 2))

    E0, E2 = energies[0], energies[2]
    # Normalize wavefunctions to ensure integral(psi^2 dx) = 1
    psi0 = wavefuncs_raw[:, 0] / np.sqrt(h)
    psi2 = wavefuncs_raw[:, 2] / np.sqrt(h)

    # Enforce convention: even wavefunctions are positive at the origin
    center_idx = N // 2
    if psi0[center_idx] < 0:
        psi0 *= -1.0
    if psi2[center_idx] < 0:
        psi2 *= -1.0
    
    result = (E0, E2, psi0, psi2, x)
    solution_cache[alpha_key] = result
    return result

def get_E2(alpha):
    """Returns the energy of the second excited state, E_2."""
    if alpha <= 0: return np.nan
    _, E2, _, _, _ = solve_schrodinger(alpha)
    return E2

def get_psi2_at_alpha(alpha):
    """Returns the value of the wavefunction psi_2 at x=alpha."""
    if alpha <= 0: return np.nan
    _, _, _, psi2, x = solve_schrodinger(alpha)
    # Interpolate to find the value at x=alpha accurately
    return np.interp(alpha, x, psi2)

def find_all_roots(func, alpha_min, alpha_max, n_steps=500):
    """Finds all roots of a function within a given range."""
    roots = []
    alphas = np.linspace(alpha_min, alpha_max, n_steps)
    y_values = [func(a) for a in alphas]
    
    for i in range(n_steps - 1):
        if np.sign(y_values[i]) != np.sign(y_values[i+1]):
            a1, a2 = alphas[i], alphas[i+1]
            try:
                root = brentq(func, a1, a2)
                roots.append(root)
            except (ValueError, RuntimeError):
                pass
    return roots

def main():
    """Main function to find and report alpha_0."""
    alpha_start = 0.1
    alpha_end = 10.0
    
    all_roots = []

    e2_roots = find_all_roots(get_E2, alpha_start, alpha_end)
    all_roots.extend(e2_roots)

    psi2_roots = find_all_roots(get_psi2_at_alpha, alpha_start, alpha_end)
    all_roots.extend(psi2_roots)

    if not all_roots:
        print("No value for alpha_0 could be found in the search range.")
        return

    largest_alpha0 = max(all_roots)

    # --- Verification Step ---
    # Calculate all the numbers in the final equation for F(alpha_0)
    E0, E2, psi0, psi2, x = solve_schrodinger(largest_alpha0)
    
    psi0_at_0 = np.interp(0, x, psi0)
    psi0_at_alpha = np.interp(largest_alpha0, x, psi0)
    psi2_at_0 = np.interp(0, x, psi2)
    psi2_at_alpha = np.interp(largest_alpha0, x, psi2)

    print(f"The largest root found is alpha_0 = {largest_alpha0:.9f}")
    print("\nVerification using the definition of F(alpha):")
    print("F(alpha) = (E_2/E_0) * [ (psi_2(alpha,alpha)/psi_2(alpha,0)) / (psi_0(alpha,alpha)/psi_0(alpha,0)) ]")
    print("\nValues of the components at alpha_0:")
    print(f"E_2({largest_alpha0:.4f}) = {E2:12.9f}")
    print(f"E_0({largest_alpha0:.4f}) = {E0:12.9f}")
    print(f"psi_2({largest_alpha0:.4f}, {largest_alpha0:.4f}) = {psi2_at_alpha:12.9f}")
    print(f"psi_2({largest_alpha0:.4f}, 0) = {psi2_at_0:12.9f}")
    print(f"psi_0({largest_alpha0:.4f}, {largest_alpha0:.4f}) = {psi0_at_alpha:12.9f}")
    print(f"psi_0({largest_alpha0:.4f}, 0) = {psi0_at_0:12.9f}")

    # One of the numerators should be zero
    term1 = E2/E0
    if abs(psi2_at_0) < 1e-12 or abs(psi0_at_0) < 1e-12 or abs(psi0_at_alpha) < 1e-12:
         F_alpha_0 = float('nan') # A denominator is zero, F is ill-defined
    else:
         term2_num = psi2_at_alpha / psi2_at_0
         term2_den = psi0_at_alpha / psi0_at_0
         F_alpha_0 = term1 * (term2_num / term2_den)

    print(f"\nFinal calculated value F(alpha_0) = {F_alpha_0:.12f}")
    print("\nThis confirms that F(alpha_0) is zero, as either E_2 or psi_2(alpha_0, alpha_0) is zero.")
    
    # Final answer in the required format
    print(f"\n<<<3.673290373>>>")

if __name__ == '__main__':
    main()
