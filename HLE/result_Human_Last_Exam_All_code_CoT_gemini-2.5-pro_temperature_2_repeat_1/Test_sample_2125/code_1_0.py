import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
import warnings

def solve_quantum_problem():
    """
    This function sets up the problem, finds the roots for F(alpha)=0,
    and prints the largest root and the corresponding equation values.
    """
    # Suppress warnings that might arise from numerical computations
    warnings.filterwarnings('ignore', category=np.RankWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning) 

    # 1. Setup Numerical Environment
    L = 20      # Half-width of the spatial grid. Needs to be large enough.
    N = 5001    # Number of grid points (odd number for a point exactly at x=0)
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    # Pre-calculate the constant kinetic energy part of the Hamiltonian matrix
    # H_kinetic = -1/2 * d^2/dx^2, discretized via finite differences.
    off_diag_T = -1.0 / (2.0 * dx**2) * np.ones(N - 1)
    diag_T = 1.0 / dx**2 * np.ones(N)
    
    # Memoization cache to store results and avoid re-computation
    memo_cache = {}

    # 2. Implement the Schrödinger Solver
    def solve_schrodinger(alpha):
        """
        Solves the TISE for a given alpha and returns the relevant
        energies (E0, E2) and wavefunctions (psi0, psi2).
        """
        if alpha in memo_cache:
            return memo_cache[alpha]

        if alpha <= 0:
            raise ValueError("alpha must be positive")

        # Potential V(x) = - (7/2)x^2 + (1/2)alpha^2*x^2 - alpha*x^4 + (1/2)x^6
        V = (-3.5 + 0.5 * alpha**2) * x**2 - alpha * x**4 + 0.5 * x**6
        
        diag_H = diag_T + V

        try:
            # Solve the eigenvalue problem, requesting the first 3 eigenstates (n=0, 1, 2)
            eigenvalues, eigenvectors = eigh_tridiagonal(
                diag_H, off_diag_T, select='i', select_range=(0, 2)
            )
        except Exception:
            return (np.nan, np.nan, None, None)

        E0, E2 = eigenvalues[0], eigenvalues[2]
        psi0_vec, psi2_vec = eigenvectors[:, 0], eigenvectors[:, 2]

        # Fix sign convention for even wavefunctions psi0 and psi2
        # to be positive at the origin for consistency.
        idx_0 = N // 2
        if psi0_vec is not None and psi0_vec[idx_0] < 0:
            psi0_vec *= -1
        if psi2_vec is not None and psi2_vec[idx_0] < 0:
            psi2_vec *= -1
            
        result = (E0, E2, psi0_vec, psi2_vec)
        memo_cache[alpha] = result
        return result

    # 3. Define the functions whose roots we need to find
    def check_E2(alpha):
        """Returns the energy of the second excited state, E2."""
        _, E2, _, _ = solve_schrodinger(alpha)
        return E2

    def check_psi2_node(alpha):
        """Returns the value of the psi_2 wavefunction at x=alpha."""
        if alpha >= L:
            return np.nan 
        _, _, _, psi2 = solve_schrodinger(alpha)
        if psi2 is None:
            return np.nan
        # Interpolate to find psi_2 at x=alpha for higher accuracy
        return np.interp(alpha, x, psi2)
        
    # 4. Find all roots in a given range
    roots = []
    def find_roots(func, a_min, a_max, steps=200):
        """Scans a range to find and store roots of a given function."""
        alphas = np.linspace(a_min, a_max, steps)
        # Vectorized evaluation of the function
        func_values = np.array([func(a) for a in alphas])
        
        for i in range(len(alphas) - 1):
            if np.isnan(func_values[i]) or np.isnan(func_values[i+1]):
                continue
            # A root exists in the interval if there's a sign change
            if np.sign(func_values[i]) != np.sign(func_values[i+1]):
                a, b = alphas[i], alphas[i+1]
                try:
                    root = brentq(func, a, b, xtol=1e-9, rtol=1e-9)
                    roots.append(root)
                except (ValueError, RuntimeError):
                    continue

    print("Solving for the largest alpha_0 where F(alpha_0) = 0.")
    print("This requires finding roots for E2(alpha) = 0 and psi2(alpha; alpha) = 0.")
    search_min, search_max = 0.1, 10.0
    print(f"Searching for roots in the range alpha = ({search_min}, {search_max})...")

    # Find roots from both conditions that make F(alpha) zero
    find_roots(check_E2, search_min, search_max)
    find_roots(check_psi2_node, search_min, search_max)
    
    # 5. Identify the largest root
    if not roots:
        print("\nCould not find any value for alpha_0 in the searched range.")
        alpha_0 = np.nan
    else:
        alpha_0 = max(roots)
        print(f"\nFound potential alpha values: {[f'{r:.5f}' for r in sorted(roots)]}")
        print(f"The largest value found is alpha_0 = {alpha_0:.5f}")

        # 6. Final output showing the numbers in the equation
        print("\nEvaluating the components of F(alpha_0) at this value:")
        
        E0, E2, psi0, psi2 = solve_schrodinger(alpha_0)
        idx_0 = N // 2
        psi0_at_0 = psi0[idx_0]
        psi2_at_0 = psi2[idx_0]
        psi0_at_alpha0 = np.interp(alpha_0, x, psi0)
        psi2_at_alpha0 = np.interp(alpha_0, x, psi2)
        
        print(f"  alpha_0 = {alpha_0:.5f}")
        print(f"  E0({alpha_0:.3f}) = {E0:.5f}")
        print(f"  E2({alpha_0:.3f}) = {E2:.5f}  <-- This term is close to zero, causing F=0.")
        print(f"  psi0({alpha_0:.3f}; x=0) = {psi0_at_0:.5f}")
        print(f"  psi0({alpha_0:.3f}; x={alpha_0:.3f}) = {psi0_at_alpha0:.5f}")
        print(f"  psi2({alpha_0:.3f}; x=0) = {psi2_at_0:.5f}")
        print(f"  psi2({alpha_0:.3f}; x={alpha_0:.3f}) = {psi2_at_alpha0:.5f}")
        
        # Check for division by zero before printing the equation
        denom_term = E0 * (psi0_at_alpha0 / psi0_at_0) * psi2_at_0
        if abs(denom_term) < 1e-9:
             print("\nFinal equation for F(alpha_0) has a denominator close to zero and is ill-defined.")
        else:
            # Reconstruct F(alpha_0)
            val = (E2 / E0) * ((psi2_at_alpha0 / psi2_at_0) / (psi0_at_alpha0 / psi0_at_0))
            print("\nFinal equation check:")
            print(f"F({alpha_0:.5f}) = ({E2:.5f} / {E0:.5f}) * (({psi2_at_alpha0:.5f} / {psi2_at_0:.5f}) / ({psi0_at_alpha0:.5f} / {psi0_at_0:.5f})) = {val:.5g}")


    # Final answer submission
    print(f"\n<<<{alpha_0:.5f}>>>")

# Execute the main function
solve_quantum_problem()