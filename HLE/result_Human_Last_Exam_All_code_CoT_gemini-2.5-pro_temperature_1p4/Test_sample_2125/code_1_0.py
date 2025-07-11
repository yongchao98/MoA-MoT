import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import interp1d
from scipy.optimize import brentq

def solve_for_alpha_zero():
    """
    This function carries out the full analysis to find the largest alpha_0
    for which F(alpha_0) = 0.
    """
    # --- Numerical Parameters ---
    # A wide grid and many points ensure accuracy, especially for larger alpha.
    X_MAX = 30.0        # Spatial grid from -X_MAX to +X_MAX
    N_POINTS = 10001    # Number of grid points

    # --- Schrödinger Equation Solver ---
    def solve_se(alpha):
        """
        Solves the 1D TISE using the finite difference method on a grid.
        Returns the first 3 energies, their corresponding wavefunctions, and the grid.
        """
        x = np.linspace(-X_MAX, X_MAX, N_POINTS)
        h = x[1] - x[0]

        # Define the potential V(x)
        V = -3.5 * x**2 + 0.5 * alpha**2 * x**2 - alpha * x**4 + 0.5 * x**6

        # Create the tridiagonal Hamiltonian matrix
        diag = 1 / h**2 + V
        off_diag = -1 / (2 * h**2) * np.ones(N_POINTS - 1)

        # Solve the eigenvalue problem. eigh_tridiagonal is efficient and returns
        # eigenvalues in ascending order.
        energies, wfs = eigh_tridiagonal(diag, off_diag)

        # Normalize the wavefunctions.
        wfs /= np.sqrt(h)

        return energies, wfs, x

    # --- Functions for Root-Finding ---
    def get_E2(alpha):
        """Returns the energy of the second excited state, E_2."""
        energies, _, _ = solve_se(alpha)
        return energies[2]

    def get_psi2_at_alpha(alpha):
        """Returns the value of the second excited state wavefunction at x = alpha."""
        if alpha >= X_MAX:
            return np.nan

        _, wfs, x = solve_se(alpha)
        psi2_vec = wfs[:, 2] # The 3rd eigenvector (index 2) is psi_2

        # The potential is symmetric, so psi_2 is an even function.
        # We enforce a consistent sign convention by making psi_2(0) > 0.
        zero_idx = N_POINTS // 2
        if psi2_vec[zero_idx] < 0:
            psi2_vec *= -1

        # Interpolate the discrete wavefunction to find its value at x = alpha.
        psi2_interp = interp1d(x, psi2_vec, kind='cubic')
        return psi2_interp(alpha)

    # --- Main Calculation ---
    all_roots = []
    
    # 1. Find root(s) for E_2(alpha) = 0
    # A scan shows E_2 crosses zero once between alpha=3 and alpha=4.
    try:
        root_E2 = brentq(get_E2, 3.0, 4.0)
        all_roots.append(root_E2)
    except ValueError:
        pass # No sign change in bracket

    # 2. Find root(s) for psi_2(alpha; alpha) = 0
    # This function oscillates, so we scan for sign changes to find all root brackets.
    alpha_scan = np.linspace(0.1, 20.0, 200)
    # Using a list comprehension and handling potential errors during the scan
    psi2_vals = []
    for a in alpha_scan:
        try:
            val = get_psi2_at_alpha(a)
            psi2_vals.append(val)
        except:
            psi2_vals.append(np.nan)

    for i in range(len(alpha_scan) - 1):
        # Check for sign change, ignoring NaNs
        if np.isfinite(psi2_vals[i]) and np.isfinite(psi2_vals[i+1]):
            if np.sign(psi2_vals[i]) != np.sign(psi2_vals[i+1]):
                a1, a2 = alpha_scan[i], alpha_scan[i+1]
                try:
                    root_psi2 = brentq(get_psi2_at_alpha, a1, a2)
                    all_roots.append(root_psi2)
                except (ValueError, RuntimeError):
                    pass

    # 3. Find the largest root and print the details
    if all_roots:
        alpha_0 = max(all_roots)

        # Calculate all components of F(alpha_0) for verification
        energies, wfs, x = solve_se(alpha_0)
        E0, E2 = energies[0], energies[2]
        
        psi0_vec, psi2_vec = wfs[:, 0], wfs[:, 2]
        zero_idx = N_POINTS // 2
        if psi0_vec[zero_idx] < 0: psi0_vec *= -1
        if psi2_vec[zero_idx] < 0: psi2_vec *= -1

        psi0_interp = interp1d(x, psi0_vec, kind='cubic')
        psi2_interp = interp1d(x, psi2_vec, kind='cubic')
        
        psi0_at_0 = psi0_interp(0)
        psi0_at_alpha0 = psi0_interp(alpha_0)
        psi2_at_0 = psi2_interp(0)
        psi2_at_alpha0 = psi2_interp(alpha_0)

        # Output the numbers in the final equation as requested
        print(f"For the largest found root alpha_0 = {alpha_0:.8f}, the components of F(alpha_0) are:")
        print(f"E_2({alpha_0:.3f}) = {E2:.8f}")
        print(f"E_0({alpha_0:.3f}) = {E0:.8f}")
        print(f"psi_2({alpha_0:.3f}; {alpha_0:.3f}) = {psi2_at_alpha0:.8e}")
        print(f"psi_2({alpha_0:.3f}; 0) = {psi2_at_0:.8f}")
        print(f"psi_0({alpha_0:.3f}; {alpha_0:.3f}) = {psi0_at_alpha0:.8f}")
        print(f"psi_0({alpha_0:.3f}; 0) = {psi0_at_0:.8f}")
        
        term1 = E2 / E0
        term2 = (psi2_at_alpha0 / psi2_at_0) / (psi0_at_alpha0 / psi0_at_0)
        F_alpha_0 = term1 * term2
        
        print(f"\nVerification: F({alpha_0:.3f}) = ({E2:.4f}/{E0:.4f}) * (({psi2_at_alpha0:.2e}/{psi2_at_0:.4f}) / ({psi0_at_alpha0:.4f}/{psi0_at_0:.4f})) = {F_alpha_0:.2e}")
        
        print(f"\nThe largest value, alpha_0, is:")
        print(f"{alpha_0}")

    else:
        print("No roots were found in the specified search range.")

if __name__ == '__main__':
    solve_for_alpha_zero()