import numpy as np
from scipy.special import lambertw
from scipy.linalg import eigh_tridiagonal

def solve_energy_difference():
    """
    Calculates the energy difference between the first and second energy levels
    for a particle in a given 3D potential well.
    """
    # --- 1. Define Constants and Parameters ---
    V0 = 15.0  # eV
    R = 3.0    # nm
    # Constant hbar^2 / (2 * m_e) in units of eV * nm^2
    hbar2_2m = 0.0380998

    # --- 2. Set up the Numerical Grid ---
    # The grid should extend far enough for the bound state wavefunctions to decay.
    r_max = 10 * R
    N = 5000  # Number of grid points for precision
    dr = r_max / N
    # The grid starts at dr to avoid the singularity at r=0.
    r_grid = np.linspace(dr, r_max, N)

    # --- 3. Define the Potential Energy Function ---
    def potential_V(r, V0_val, R_val):
        """Calculates the potential V(r) on the given grid."""
        V = np.zeros_like(r)
        
        # Define masks for the two regions of the potential
        mask_lt_R = r < R_val
        mask_ge_R = ~mask_lt_R
        
        r_lt_R = r[mask_lt_R]
        r_ge_R = r[mask_ge_R]
        
        # For 0 <= r < R (repulsive core)
        if r_lt_R.size > 0:
            # The argument to lambertw is always real and positive.
            # np.real is used to ensure the output is float.
            V[mask_lt_R] = np.sqrt(V0_val + np.real(lambertw(np.exp(r_lt_R - R_val))))
        
        # For r >= R (potential well)
        if r_ge_R.size > 0:
            V[mask_ge_R] = np.sqrt(V0_val * (1.0 - (R_val / r_ge_R)**2))
            
        return V

    # --- 4. Solve the Schrödinger Equation Numerically ---
    def get_energies(l_val, r, dr_val, V_r):
        """
        Solves the radial Schrödinger equation for a given angular momentum l.
        """
        # Define the effective potential, including the centrifugal barrier
        V_eff = V_r + l_val * (l_val + 1) * hbar2_2m / r**2
        
        # Set up the Hamiltonian matrix components
        diag = 2 * hbar2_2m / dr_val**2 + V_eff
        off_diag = -hbar2_2m / dr_val**2 * np.ones(N - 1)
        
        # Solve the eigenvalue problem for the tridiagonal matrix
        # This returns eigenvalues (energies) sorted in ascending order.
        eigenvalues, _ = eigh_tridiagonal(diag, off_diag)
        
        # Filter for bound states where E < V(infinity)
        V_inf = np.sqrt(V0)
        bound_eigenvalues = eigenvalues[eigenvalues < V_inf]
        
        return bound_eigenvalues

    # Calculate the potential V(r) on the grid
    V_r = potential_V(r_grid, V0, R)
    
    # Calculate energy levels for l=0 (s-states) and l=1 (p-states)
    energies_l0 = get_energies(0, r_grid, dr, V_r)
    energies_l1 = get_energies(1, r_grid, dr, V_r)

    # --- 5. Identify E1 and E2 and Calculate the Difference ---
    # The ground state (E1) is the lowest energy state, which is the 1s state.
    E1 = energies_l0[0]
    
    # The first excited state (E2) is the second lowest energy state overall.
    # We must compare the second s-state (2s) and the first p-state (1p).
    E_2s = energies_l0[1]
    E_1p = energies_l1[0]
    E2 = min(E_2s, E_1p)
    
    # Calculate the energy difference
    delta_E = E2 - E1

    # --- 6. Print the Final Result ---
    print("The energy difference ΔE = E2 - E1 is calculated as follows:")
    print(f"Ground state energy E1 = {E1:.4f} eV")
    print(f"First excited state energy E2 = {E2:.4f} eV")
    print(f"ΔE = {E2:.4f} eV - {E1:.4f} eV = {delta_E:.4f} eV")
    
    # Final answer in the required format
    print(f"\n<<<{delta_E:.4f}>>>")

# Execute the function to get the answer
solve_energy_difference()