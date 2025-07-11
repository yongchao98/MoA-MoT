import numpy as np
from scipy.special import lambertw
from scipy.constants import hbar, m_e, e

def solve_quantum_well_energies():
    """
    Calculates the energy difference between the first two energy levels of a particle
    in a 3D potential well by numerically solving the Schrödinger equation.
    """
    # Step 1: Define constants and parameters
    V0_eV = 15.0  # Potential energy at infinity in eV
    R_nm = 3.0    # Well radius in nm
    
    # Convert parameters to SI units for calculation
    V0_si = V0_eV * e
    R_si = R_nm * 1e-9
    m = m_e

    # Step 2: Set up the numerical grid for the radial coordinate r
    # A fine grid extending well beyond the well radius is used for accuracy.
    num_points = 4000
    r_max = 10 * R_si  # Solve up to 30 nm
    r, dr = np.linspace(0, r_max, num_points, retstep=True)
    r[0] = 1e-12 # Start grid very close to zero to avoid any division-by-zero errors.
    
    # Step 3: Define and compute the base potential V(r) on the grid
    # This uses the physically corrected version of the potential.
    def get_potential(r_coords, V0, R):
        x = r_coords / R
        V = np.zeros_like(r_coords)
        
        mask_lt_R = (x < 1) & (x > 0)
        mask_ge_R = (x >= 1)
        
        # For r < R
        W_arg = np.exp(x[mask_lt_R] - 1)
        W_val = lambertw(W_arg).real
        V[mask_lt_R] = V0 * np.sqrt(1 + W_val)
        
        # For r >= R
        V[mask_ge_R] = V0 * np.sqrt(1 - (x[mask_ge_R])**-2)
        
        return V

    V_base_grid = get_potential(r, V0_si, R_si)
    
    # Pre-calculate constant terms for the Hamiltonian matrix
    T_const = hbar**2 / (m * dr**2)
    U_const = -hbar**2 / (2 * m * dr**2)
    off_diagonal_elements = np.full(num_points - 1, U_const)
    
    all_energies_eV = []
    
    # Step 4: Solve for energy levels for different angular momentum l
    # We check l=0, 1, 2 to find the lowest energy states of the system.
    for l in range(3):
        # Add centrifugal potential term for l > 0
        V_eff = V_base_grid
        if l > 0:
            centrifugal_term = hbar**2 * l * (l + 1) / (2 * m * r**2)
            V_eff = V_base_grid + centrifugal_term

        # Construct Hamiltonian matrix for this l
        diagonal_elements = T_const + V_eff
        H = np.diag(diagonal_elements) + np.diag(off_diagonal_elements, k=1) + np.diag(off_diagonal_elements, k=-1)
        
        # Find eigenvalues (energies)
        eigenvalues_J = np.linalg.eigh(H)[0]
        
        # Convert to eV and filter for bound states (E < V0)
        eigenvalues_eV = eigenvalues_J / e
        bound_states_eV = eigenvalues_eV[eigenvalues_eV < V0_eV]
        
        # Store the found energy levels
        all_energies_eV.extend(bound_states_eV)

    # Step 5: Determine the first and second energy levels
    # Sort all unique energies found across different l values
    sorted_unique_energies = sorted(list(set(all_energies_eV)))
    
    if len(sorted_unique_energies) < 2:
        print("Could not find at least two distinct bound energy levels.")
        return

    E1_eV = sorted_unique_energies[0]
    E2_eV = sorted_unique_energies[1]
    
    # Step 6: Calculate the energy difference and print the result
    delta_E_eV = E2_eV - E1_eV
    
    print(f"The first energy level (E1) is: {E1_eV:.4f} eV")
    print(f"The second energy level (E2) is: {E2_eV:.4f} eV")
    print(f"The energy difference is Delta E = E2 - E1 = {E2_eV:.4f} - {E1_eV:.4f} = {delta_E_eV:.4f} eV")
    
    return delta_E_eV

# Execute the main function
final_answer = solve_quantum_well_energies()
print(f"\n<<<{final_answer:.4f}>>>")