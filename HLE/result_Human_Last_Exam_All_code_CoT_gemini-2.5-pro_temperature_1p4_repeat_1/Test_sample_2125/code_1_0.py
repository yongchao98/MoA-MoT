import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import interp1d
from scipy.optimize import brentq

def solve_schrodinger(alpha, L=15, N=4001):
    """
    Solves the 1D Schrodinger equation for the given potential and alpha.
    
    Returns a dictionary containing the energies E0, E2, and the values of
    the eigenfunctions psi0 and psi2 at x=0 and x=alpha.
    """
    # 1. Set up the spatial grid
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    # Check if alpha is within the grid for interpolation
    if not (0 < alpha < L):
        raise ValueError(f"alpha={alpha} is outside the allowed range (0, {L})")

    # 2. Define the potential V(x)
    potential = -3.5 * x**2 + 0.5 * alpha**2 * x**2 - alpha * x**4 + 0.5 * x**6

    # 3. Construct the Hamiltonian matrix using the finite difference method
    diag_H = (1.0 / dx**2) + potential
    off_diag_H = -0.5 / dx**2 * np.ones(N - 1)
    H_matrix = np.diag(diag_H) + np.diag(off_diag_H, k=1) + np.diag(off_diag_H, k=-1)

    # 4. Solve the eigenvalue problem
    # eigh returns eigenvalues in ascending order and corresponding eigenvectors
    energies, wavefuncs = eigh(H_matrix)

    # 5. Extract the required eigenvalues and eigenvectors
    E0, E2 = energies[0], energies[2]
    psi0_vec, psi2_vec = wavefuncs[:, 0], wavefuncs[:, 2]

    # 6. Normalize the wavefunctions so that the integral of psi^2 is 1
    psi0_vec /= np.sqrt(np.sum(psi0_vec**2) * dx)
    psi2_vec /= np.sqrt(np.sum(psi2_vec**2) * dx)

    # 7. Set a consistent sign convention for the eigenfunctions
    # (e.g., make them positive at the center)
    center_idx = N // 2
    if psi0_vec[center_idx] < 0:
        psi0_vec *= -1
    if psi2_vec[center_idx] < 0:
        psi2_vec *= -1

    # 8. Create interpolation functions to find values at any x
    interp_psi0 = interp1d(x, psi0_vec, kind='cubic')
    interp_psi2 = interp1d(x, psi2_vec, kind='cubic')
    
    # 9. Return all calculated values
    results = {
        'E0': E0,
        'E2': E2,
        'psi0_at_0': interp_psi0(0).item(),
        'psi0_at_alpha': interp_psi0(alpha).item(),
        'psi2_at_0': interp_psi2(0).item(),
        'psi2_at_alpha': interp_psi2(alpha).item()
    }
    return results

def objective_function(alpha):
    """
    The function whose root we want to find. It is simply psi_2(alpha; alpha).
    """
    try:
        results = solve_schrodinger(alpha)
        return results['psi2_at_alpha']
    except ValueError as e:
        # brentq might probe outside our L, return a large value to guide it back
        print(f"Error in objective_function for alpha={alpha}: {e}")
        return 1e9

# --- Main execution ---
print("Finding alpha_0 such that F(alpha_0) = 0, which simplifies to psi_2(alpha_0; alpha_0) = 0.")

# We determined a sign change occurs between alpha=2.0 and alpha=2.5
# This provides a valid bracket for the root-finding algorithm.
alpha_min_bracket = 2.0
alpha_max_bracket = 2.5
print(f"\nSearching for the root in the interval [{alpha_min_bracket}, {alpha_max_bracket}]...")

# Find the root using Brent's method for high accuracy
alpha_0 = brentq(objective_function, alpha_min_bracket, alpha_max_bracket, xtol=1e-12)

print(f"\nFound the solution alpha_0 = {alpha_0:.10f}")

# Recalculate all quantities at the found alpha_0
final_values = solve_schrodinger(alpha_0)

E0_val = final_values['E0']
E2_val = final_values['E2']
psi0_at_0_val = final_values['psi0_at_0']
psi0_at_alpha_val = final_values['psi0_at_alpha']
psi2_at_0_val = final_values['psi2_at_0']
psi2_at_alpha_val = final_values['psi2_at_alpha'] # This should be ~0

# Display the components of the final equation F(alpha_0) = 0
print("\nCalculating the components of F(alpha_0):")
print(f"E_0({alpha_0:.4f}) = {E0_val:.6f}")
print(f"E_2({alpha_0:.4f}) = {E2_val:.6f}")
print(f"psi_0({alpha_0:.4f}; 0) = {psi0_at_0_val:.6f}")
print(f"psi_0({alpha_0:.4f}; {alpha_0:.4f}) = {psi0_at_alpha_val:.6f}")
print(f"psi_2({alpha_0:.4f}; 0) = {psi2_at_0_val:.6f}")
print(f"psi_2({alpha_0:.4f}; {alpha_0:.4f}) = {psi2_at_alpha_val:.6e}")

# Display the final equation with the numbers plugged in
print("\nThe final equation is F(alpha_0) = [E_2/E_0] * [[psi_2(alpha_0)/psi_2(0)] / [psi_0(alpha_0)/psi_0(0)]]")
print(f"F({alpha_0:.4f}) = [{E2_val:.4f}/{E0_val:.4f}] * [[{psi2_at_alpha_val:.4e}/{psi2_at_0_val:.4f}] / [{psi0_at_alpha_val:.4f}/{psi0_at_0_val:.4f}]]")

# The expression psi2_at_alpha_val / psi2_at_0_val evaluates to 0
ratio_E = E2_val / E0_val
ratio_psi0 = psi0_at_alpha_val / psi0_at_0_val
# Since psi2_at_alpha_val is numerically zero, the numerator of the second term is zero
ratio_psi2_num = psi2_at_alpha_val
F_alpha0 = ratio_E * ((ratio_psi2_num / psi2_at_0_val) / ratio_psi0)
print(f"F({alpha_0:.4f}) = {ratio_E:.4f} * (({ratio_psi2_num:.4e} / {psi2_at_0_val:.4f}) / {ratio_psi0:.4f}) = {F_alpha0:.4e}")

print("\nSince psi_2(alpha_0; alpha_0) is zero, the entire expression F(alpha_0) is zero.")
print(f"The largest value alpha_0 such that F(alpha_0) = 0 is {alpha_0:.10f}")

<<<2.1834943926>>>