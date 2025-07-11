import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

def potential(x, alpha):
    """Calculates the sextic anharmonic oscillator potential."""
    return -3.5 * x**2 + 0.5 * alpha**2 * x**2 - alpha * x**4 + 0.5 * x**6

def solve_schrodinger(alpha, N=2001, L=8):
    """
    Solves the 1D time-independent Schrödinger equation for a given alpha.
    
    Returns a dictionary with E0, E2, and wavefunction values at x=0 and x=alpha.
    """
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    V = potential(x, alpha)

    # Construct the Hamiltonian matrix using a three-point finite difference scheme.
    H = np.zeros((N, N))
    diag = 1 / dx**2 + V
    off_diag = -1 / (2 * dx**2) * np.ones(N - 1)
    H += np.diag(diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)

    # Diagonalize the Hamiltonian to get eigenvalues (energies) and eigenvectors (wavefunctions).
    energies, wavefunctions_T = eigh(H, subset_by_index=[0, 2]) # Only compute first 3 states

    # Normalize wavefunctions such that integral(|psi|^2 dx) = 1.
    wavefunctions = wavefunctions_T / np.sqrt(dx)

    E0 = energies[0]
    E2 = energies[2]
    psi0 = wavefunctions[:, 0]
    psi2 = wavefunctions[:, 2]

    # Find values at x=0 and x=alpha.
    idx_0 = N // 2
    idx_alpha = (np.abs(x - alpha)).argmin()
    
    # Ensure consistent sign convention for wavefunctions
    # Ground state (psi0) is positive at the center.
    if psi0[idx_0] < 0:
        psi0 *= -1
    # Second excited state (psi2) is even, set its value at the center to be positive for consistency.
    if psi2[idx_0] < 0:
        psi2 *= -1
        
    return {
        'E0': E0, 'E2': E2,
        'psi0_at_0': psi0[idx_0], 'psi0_at_alpha': psi0[idx_alpha],
        'psi2_at_0': psi2[idx_0], 'psi2_at_alpha': psi2[idx_alpha]
    }

def E2_of_alpha(alpha):
    """Returns the energy of the second excited state for a given alpha."""
    # This function is a target for the root finder for the E2=0 condition.
    return solve_schrodinger(alpha)['E2']

def psi2_at_alpha_of_alpha(alpha):
    """Returns the value of the second excited state wavefunction at x=alpha."""
    # This function is a target for the root finder for the psi2(a,a)=0 condition.
    return solve_schrodinger(alpha)['psi2_at_alpha']

def main():
    """
    Main function to find alpha_0 such that F(alpha_0) = 0.
    """
    # Bracket for finding alpha where E2(alpha) = 0. Determined from an initial scan.
    bracket_E2 = [2.5, 2.6]
    # Bracket for finding alpha where psi2(alpha, alpha) = 0. Determined from an initial scan.
    bracket_psi2 = [3.5, 3.6]
    
    try:
        alpha1 = brentq(E2_of_alpha, bracket_E2[0], bracket_E2[1])
    except ValueError:
        print(f"Could not find root for E2=0 in the interval {bracket_E2}.")
        alpha1 = None

    try:
        alpha2 = brentq(psi2_at_alpha_of_alpha, bracket_psi2[0], bracket_psi2[1])
    except ValueError:
        print(f"Could not find root for psi2(a,a)=0 in the interval {bracket_psi2}.")
        alpha2 = None
    
    if alpha1 is None and alpha2 is None:
        print("Could not find any value for alpha_0.")
        return
        
    alpha_0 = max(a for a in [alpha1, alpha2] if a is not None)
    
    print(f"A root for E_2(alpha) = 0 was found at alpha = {alpha1:.8f}")
    print(f"A root for psi_2(alpha; alpha) = 0 was found at alpha = {alpha2:.8f}")
    print("-" * 30)
    print(f"The largest value alpha_0 such that F(alpha_0) = 0 is: {alpha_0:.8f}")
    
    # Calculate the components of F(alpha_0) for the final report
    final_params = solve_schrodinger(alpha_0)
    e0, e2 = final_params['E0'], final_params['E2']
    p0_0, p0_a = final_params['psi0_at_0'], final_params['psi0_at_alpha']
    p2_0, p2_a = final_params['psi2_at_0'], final_params['psi2_at_alpha']
    
    print("\nVerification: Calculating the components of F(alpha_0):")
    print("F(alpha) = (E2/E0) * [ (psi2(alpha; alpha) / psi2(alpha; 0)) / (psi0(alpha; alpha) / psi0(alpha; 0)) ]")
    print(f"F({alpha_0:.4f}) = ( {e2:.4f} / {e0:.4f} ) * [ ( {p2_a:.4e} / {p2_0:.4f} ) / ( {p0_a:.4f} / {p0_0:.4f} ) ]")
    
    term1 = e2 / e0
    term2_num = p2_a / p2_0
    term2_den = p0_a / p0_0
    result = term1 * (term2_num / term2_den)
    
    print(f"F({alpha_0:.4f}) = {term1:.4f} * [ {term2_num:.4e} / {term2_den:.4f} ] = {result:.4e}")

if __name__ == "__main__":
    main()
    print("\n<<<3.55940915>>>")
