import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

def solve_schrodinger(alpha, x_grid):
    """
    Solves the 1D time-independent Schrodinger equation for the given potential.

    Args:
        alpha (float): The parameter in the potential V(x).
        x_grid (np.array): The spatial grid for discretization.

    Returns:
        tuple: A tuple containing (E0, E2, psi0, psi2).
               E0, E2: Ground state and 2nd excited state energies.
               psi0, psi2: Corresponding normalized wavefunctions on the grid.
    """
    N = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    
    # Potential V(x)
    V = (0.5 * alpha**2 - 3.5) * x_grid**2 - alpha * x_grid**4 + 0.5 * x_grid**6
    
    # Kinetic energy operator (finite difference)
    T_diag = np.ones(N) / (dx**2)
    T_off_diag = -0.5 * np.ones(N - 1) / (dx**2)
    
    # Hamiltonian matrix
    H = np.diag(T_diag + V) + np.diag(T_off_diag, k=1) + np.diag(T_off_diag, k=-1)
    
    # Solve for the 3 lowest eigenvalues and eigenvectors
    # eigh is efficient for symmetric matrices
    eigenvalues, eigenvectors = eigh(H, subset_by_index=[0, 2])
    
    E0, E2 = eigenvalues[0], eigenvalues[2]
    psi0, psi2 = eigenvectors[:, 0], eigenvectors[:, 2]

    # Normalize wavefunctions such that integral(psi^2)dx = 1
    psi0 /= np.sqrt(np.sum(psi0**2) * dx)
    psi2 /= np.sqrt(np.sum(psi2**2) * dx)
    
    # Enforce consistent phase (positive at x=0 for these even states)
    mid_idx = N // 2
    if psi0[mid_idx] < 0:
        psi0 *= -1
    if psi2[mid_idx] < 0:
        psi2 *= -1
        
    return E0, E2, psi0, psi2

# Cached solver to avoid redundant computations during root finding
_cache = {}
def cached_solver(alpha, x_grid):
    alpha_key = round(alpha, 12) # Use rounded alpha as key to handle precision issues
    if alpha_key in _cache:
        return _cache[alpha_key]
    result = solve_schrodinger(alpha, x_grid)
    _cache[alpha_key] = result
    return result

# Define functions for root finding
def get_E2(alpha, x_grid):
    _, E2, _, _ = cached_solver(alpha, x_grid)
    return E2

def get_psi2_at_alpha(alpha, x_grid):
    _, _, _, psi2 = cached_solver(alpha, x_grid)
    spline = CubicSpline(x_grid, psi2)
    return float(spline(alpha))

def main():
    """Main function to execute the plan and find alpha_0."""
    # Setup the numerical grid
    N = 2001
    L = 12.0
    x_grid = np.linspace(-L, L, N)

    # --- Find root for E2(alpha) = 0 ---
    # Based on analysis, a root is expected between alpha=2.3 and 2.4
    a_e2, b_e2 = 2.3, 2.4
    print(f"Searching for root of E_2(alpha) = 0 in interval [{a_e2}, {b_e2}]...")
    root_E2 = brentq(get_E2, a_e2, b_e2, args=(x_grid,))
    print(f"Found root for E_2(alpha) = 0 at alpha = {root_E2:.8f}\n")

    # --- Find root for psi2(alpha; alpha) = 0 ---
    # Based on analysis, a root is expected between alpha=1.6 and 1.7
    a_psi, b_psi = 1.6, 1.7
    print(f"Searching for root of psi_2(alpha; alpha) = 0 in interval [{a_psi}, {b_psi}]...")
    root_psi2 = brentq(get_psi2_at_alpha, a_psi, b_psi, args=(x_grid,))
    print(f"Found root for psi_2(alpha; alpha) = 0 at alpha = {root_psi2:.8f}\n")
    
    # --- Determine the largest root, alpha_0 ---
    alpha_0 = max(root_E2, root_psi2)
    print(f"The candidate values for alpha_0 are {root_E2:.8f} and {root_psi2:.8f}.")
    print(f"The largest value is alpha_0 = {alpha_0:.8f}\n")

    # --- Verify the solution by calculating F(alpha_0) ---
    print("Verifying the solution by calculating the components of F(alpha_0):")
    E0_val, E2_val, psi0_vec, psi2_vec = cached_solver(alpha_0, x_grid)
    
    # Interpolate wavefunctions at required points
    spline_psi0 = CubicSpline(x_grid, psi0_vec)
    spline_psi2 = CubicSpline(x_grid, psi2_vec)
    
    psi0_at_0 = float(spline_psi0(0))
    psi0_at_a0 = float(spline_psi0(alpha_0))
    psi2_at_0 = float(spline_psi2(0))
    psi2_at_a0 = float(spline_psi2(alpha_0))
    
    # Assemble F(alpha_0)
    term1 = E2_val / E0_val
    term2_num = psi2_at_a0 / psi2_at_0
    term2_den = psi0_at_a0 / psi0_at_0
    
    # Handle the case where the numerator of a fraction is zero
    if abs(psi2_at_a0) < 1e-9: # This term makes F(alpha) zero
        term2 = 0
    else:
        term2 = term2_num / term2_den
    
    F_alpha0 = term1 * term2
    
    # Print the final equation with numbers
    print(f"For alpha_0 = {alpha_0:.8f}:")
    print(f"  E_0({alpha_0:.4f}) = {E0_val:.8f}")
    print(f"  E_2({alpha_0:.4f}) = {E2_val:.8f}")
    print(f"  psi_0({alpha_0:.4f}; 0) = {psi0_at_0:.8f}")
    print(f"  psi_0({alpha_0:.4f}; {alpha_0:.4f}) = {psi0_at_a0:.8f}")
    print(f"  psi_2({alpha_0:.4f}; 0) = {psi2_at_0:.8f}")
    print(f"  psi_2({alpha_0:.4f}; {alpha_0:.4f}) = {psi2_at_a0:.8f}\n")
    
    print("Final equation F(alpha_0):")
    print(f"F({alpha_0:.4f}) = ({E2_val:.8f} / {E0_val:.8f}) * (({psi2_at_a0:.8f} / {psi2_at_0:.8f}) / ({psi0_at_a0:.8f} / {psi0_at_0:.8f}))")
    print(f"F({alpha_0:.4f}) = {F_alpha0:.8e}")
    print("\nThe largest value, alpha_0, such that F(alpha_0) = 0 is approximately:")
    print(f"{alpha_0:.8f}")

if __name__ == "__main__":
    main()
    # Outputting the final answer in the requested format.
    # Recalculating here to be sure, although main() prints it.
    N_final = 2001
    L_final = 12.0
    x_grid_final = np.linspace(-L_final, L_final, N_final)
    root_E2_final = brentq(get_E2, 2.3, 2.4, args=(x_grid_final,))
    root_psi2_final = brentq(get_psi2_at_alpha, 1.6, 1.7, args=(x_grid_final,))
    alpha_0_final = max(root_E2_final, root_psi2_final)
    print(f"\n<<<{alpha_0_final:.8f}>>>")
