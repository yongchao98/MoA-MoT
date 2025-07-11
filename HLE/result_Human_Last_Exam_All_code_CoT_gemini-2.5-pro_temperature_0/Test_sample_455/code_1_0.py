import numpy as np
from scipy.special import lambertw
from scipy.linalg import eigh_tridiagonal

def solve_energy_levels():
    """
    Calculates the energy difference between the first and second energy levels
    for a particle in a specified 3D potential well.
    """
    # 1. Define constants and parameters in consistent units (eV, nm)
    V0 = 15.0  # Potential parameter in eV
    R = 3.0    # Radius of the well in nm
    m_particle_kg = 9.11e-31 # Mass of the particle in kg

    # Use electron properties for unit conversion
    mc2_electron = 510998.95 # Electron rest mass energy in eV
    hbar_c = 197.327         # h-bar * c in eV*nm

    # Calculate the kinetic energy coefficient hbar^2 / (2m) in eV*nm^2
    hbar2_over_2m = hbar_c**2 / (2 * mc2_electron)

    # 2. Set up the numerical grid for the radial coordinate r
    # A higher number of points (N) and a sufficiently large box (r_max) ensure accuracy.
    N = 4000
    r_max = 8 * R  # Maximum radius for the simulation box
    
    # Create a grid for r from a small value > 0 up to r_max.
    # The radial function u(r) is zero at r=0, so we can exclude this point.
    r_grid, dr = np.linspace(0, r_max, N, retstep=True)
    r_grid = r_grid[1:] # Grid points from dr to r_max

    # 3. Define and calculate the potential energy V(r) on the grid
    def potential_V(r_vals, V0_val, R_val):
        """Calculates the potential V(r) based on the given piecewise function for V^2(r)."""
        V_squared = np.zeros_like(r_vals)

        # Create boolean masks for the two conditions on r
        mask_lt_R = r_vals < R_val
        mask_ge_R = r_vals >= R_val

        # Case 1: 0 <= r < R
        # np.real is used as lambertw can return complex numbers for some inputs
        V_squared[mask_lt_R] = V0_val + np.real(lambertw(np.exp(r_vals[mask_lt_R] - R_val)))

        # Case 2: r >= R
        V_squared[mask_ge_R] = V0_val * (1 - (r_vals[mask_ge_R] / R_val)**(-2))

        # V(r) is the square root of V^2(r). Ensure non-negative before sqrt.
        V_squared[V_squared < 0] = 0
        return np.sqrt(V_squared)

    V_on_grid = potential_V(r_grid, V0, R)

    # 4. Construct the tridiagonal Hamiltonian matrix
    # The kinetic term from the finite difference approximation of -d^2/dr^2
    kinetic_term = hbar2_over_2m / dr**2

    # Main diagonal of the Hamiltonian matrix: 2*T + V
    main_diag = 2 * kinetic_term + V_on_grid

    # Off-diagonal elements of the Hamiltonian matrix: -T
    off_diag = -kinetic_term * np.ones(len(r_grid) - 1)

    # 5. Solve the eigenvalue problem for the two lowest energies
    # eigh_tridiagonal is efficient for finding eigenvalues of a tridiagonal matrix.
    # We request the first 2 eigenvalues (indices 0 and 1).
    eigenvalues = eigh_tridiagonal(main_diag, off_diag, select='i', select_range=(0, 1))[0]

    # The eigenvalues are the energy levels E1 and E2
    E1 = eigenvalues[0]
    E2 = eigenvalues[1]

    # 6. Calculate the energy difference and print the final equation
    delta_E = E2 - E1

    print("The first energy level is E1 = {:.4f} eV".format(E1))
    print("The second energy level is E2 = {:.4f} eV".format(E2))
    print("\nThe energy difference is calculated as:")
    print("ΔE = E2 - E1 = {:.4f} eV - {:.4f} eV = {:.4f} eV".format(E2, E1, delta_E))
    
    # Return the final numerical answer for the submission format
    return delta_E

# Execute the function and print the final answer in the required format
final_answer = solve_energy_levels()
print(f"\n<<<{final_answer:.4f}>>>")
