import numpy as np
from scipy.special import lambertw
from scipy.linalg import eigh_tridiagonal

def solve_energy_levels():
    """
    This function calculates the energy difference between the first two energy
    levels of a quantum particle in a given 3D potential well.

    The method involves numerically solving the time-independent radial
    Schrödinger equation for s-wave states (l=0).
    """

    # --- Step 1: Define Parameters and Constants ---
    # The problem provides the following parameters:
    V0_eV = 15.0  # Potential parameter in eV
    R_nm = 3.0    # Well radius in nm
    m = 9.11e-31  # Mass of the particle (electron) in kg

    # Physical constants in SI units
    hbar = 1.054571817e-34  # Reduced Planck constant in J*s
    e = 1.602176634e-19    # Elementary charge in C

    # Convert parameters to SI units for calculations
    R_si = R_nm * 1e-9  # Radius in meters

    # --- Step 2: Interpret and Implement the Potential Function ---
    # The potential V(r) is derived from the given V^2(r). We take the
    # positive square root. The expression V^2(r) is assumed to be in units
    # of eV^2. The term (r - R) in the exponent is interpreted as using the
    # numerical values of r and R in nanometers to ensure the argument is
    # dimensionless.

    def get_potential_V_joules(r_si_grid, R_si_val, V0_eV_val):
        """
        Calculates the potential energy V(r) in Joules on a grid of r values.
        """
        # Convert SI radius grid to nanometers for use in the formula
        r_nm_grid = r_si_grid * 1e9
        R_nm_val = R_si_val * 1e9

        # Initialize array for the potential V^2 in eV^2
        V2_ev2 = np.zeros_like(r_si_grid)

        # Create boolean masks for the two regions of the potential
        mask_inside = r_si_grid < R_si_val
        mask_outside = ~mask_inside

        # Region 1: 0 <= r < R
        # V^2(r) [eV^2] = V0 + W(e^(r[nm] - R[nm]))
        r_nm_inside = r_nm_grid[mask_inside]
        arg_W = np.exp(r_nm_inside - R_nm_val)
        # The Lambert W function can be complex, but for a positive argument,
        # the principal branch (k=0) is real.
        V2_ev2[mask_inside] = V0_eV_val + np.real(lambertw(arg_W))

        # Region 2: r >= R
        # V^2(r) [eV^2] = V0 * (1 - (r/R)^-2) = V0 * (1 - (R/r)^2)
        r_si_outside = r_si_grid[mask_outside]
        V2_ev2[mask_outside] = V0_eV_val * (1 - (R_si_val / r_si_outside)**2)
        
        # Calculate V(r) in eV by taking the positive square root
        # Handle potential negative values from numerical precision, though unlikely here
        V_ev = np.sqrt(np.maximum(0, V2_ev2))

        # Convert potential from eV to Joules
        return V_ev * e

    # --- Step 3: Set up Numerical Grid and Hamiltonian ---
    # We solve the radial Schrödinger equation for u(r) = r*R(r) and l=0:
    # (-ħ²/2m) * u''(r) + V(r)*u(r) = E*u(r)
    # This is discretized on a grid and becomes a matrix eigenvalue problem.

    # Set up the spatial grid for r
    N = 4000  # Number of grid points (for accuracy)
    r_max = 8 * R_si  # Grid extends well beyond the potential well radius
    dr = r_max / (N + 1)  # Grid spacing
    r_grid = np.linspace(dr, r_max, N) # Grid from dr to r_max

    # Calculate the potential on this grid
    V_J_grid = get_potential_V_joules(r_grid, R_si, V0_eV)

    # Construct the tridiagonal Hamiltonian matrix
    # The kinetic part comes from the finite-difference approximation of u''(r)
    kinetic_off_diag = -hbar**2 / (2 * m * dr**2)
    kinetic_diag = -2 * kinetic_off_diag

    # The full Hamiltonian matrix H = T + V
    H_diag = kinetic_diag + V_J_grid
    H_offdiag = np.full(N - 1, kinetic_off_diag)

    # --- Step 4: Solve for Eigenvalues ---
    # Use scipy's efficient tridiagonal solver. It returns sorted eigenvalues.
    eigenvalues_J, _ = eigh_tridiagonal(H_diag, H_offdiag)

    # The two lowest energy levels are the first two eigenvalues
    E1_J = eigenvalues_J[0]
    E2_J = eigenvalues_J[1]

    # --- Step 5: Calculate and Print the Result ---
    # Convert energies from Joules to eV
    E1_eV = E1_J / e
    E2_eV = E2_J / e
    delta_E_eV = E2_eV - E1_eV

    print("Calculation of the energy difference ΔE = E2 - E1:")
    print("-" * 50)
    print(f"The first energy level (ground state) is E1 = {E1_eV:.4f} eV")
    print(f"The second energy level (first excited state) is E2 = {E2_eV:.4f} eV")
    print("-" * 50)
    print("The energy difference is:")
    print(f"ΔE = E2 - E1 = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")
    
    # Return the final numerical answer in the required format
    return f"<<<{delta_E_eV:.4f}>>>"

# Run the calculation and print the final answer
final_answer = solve_energy_levels()
print(final_answer)