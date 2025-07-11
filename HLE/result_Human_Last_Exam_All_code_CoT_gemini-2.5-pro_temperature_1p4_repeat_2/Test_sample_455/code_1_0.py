import numpy as np
from scipy.special import lambertw
from scipy.linalg import eigh_tridiagonal

def solve_energy_levels():
    """
    Solves for the energy levels of a quantum particle in a specified potential well.
    """
    # Step 1: Define Constants and Parameters in mixed units (eV, nm)
    V0_val = 15.0  # Numerical value for V0, used in the potential formula
    R_nm = 3.0     # Radius of the well in nm

    # The kinetic energy constant ħ²/(2m) is needed in units of eV·nm².
    # (ħ^2 / (2 * m_e)) = 0.0380998 eV·nm^2
    hbar_sq_over_2m_eV_nm2 = 0.0380998

    # Step 2: Define the Potential Energy Function V(r)
    def potential_in_eV(r_in_nm, V0, R):
        """
        Calculates the potential V(r) in units of eV for a given radius r in nm.
        """
        r = np.asarray(r_in_nm, dtype=float)
        V_squared = np.zeros_like(r)

        # Region 1: 0 <= r < R
        mask_lt_R = r < R
        r_lt_R = r[mask_lt_R]
        arg_W = np.exp(r_lt_R - R)
        V_squared[mask_lt_R] = V0 + np.real(lambertw(arg_W))

        # Region 2: r >= R
        mask_ge_R = r >= R
        r_ge_R = r[mask_ge_R]
        V_squared[mask_ge_R] = V0 * (1 - (R / r_ge_R)**2)

        # Return V(r) = sqrt(V^2(r)). V_squared is non-negative for all r >= 0.
        return np.sqrt(V_squared)

    # Step 3: Set up the Numerical Grid
    # A fine grid is needed to accurately model the sharp potential well at r=R.
    num_points = 4000
    r_max_nm = 5.0 * R_nm  # Grid extends to 5R to ensure wavefunction decay

    # Create the radial grid from a small value up to r_max.
    # The wavefunction u(r) = rR(r) is zero at r=0, so we exclude this point.
    r_grid, dr_nm = np.linspace(0, r_max_nm, num_points + 1, retstep=True)
    r_grid = r_grid[1:]  # Grid points from dr to r_max

    # Step 4: Construct the Tridiagonal Hamiltonian Matrix H
    # Using the finite difference approximation, H has diagonal and off-diagonal elements.
    # Kinetic term T = ħ² / (2 * m * dr²)
    kinetic_term_T_eV = hbar_sq_over_2m_eV_nm2 / (dr_nm**2)

    # Calculate potential V(r) on the grid
    potential_grid_eV = potential_in_eV(r_grid, V0_val, R_nm)

    # Diagonal elements of H: H_ii = 2*T + V_i
    main_diagonal = 2 * kinetic_term_T_eV + potential_grid_eV

    # Off-diagonal elements of H: H_{i, i+1} = -T
    off_diagonal = -kinetic_term_T_eV * np.ones(num_points - 1)

    # Step 5: Solve the Eigenvalue Problem
    # eigh_tridiagonal efficiently finds eigenvalues for a symmetric tridiagonal matrix.
    # The eigenvalues are the energy levels, sorted in ascending order.
    eigenvalues_eV, _ = eigh_tridiagonal(main_diagonal, off_diagonal)

    # Step 6: Extract the First Two Energy Levels
    E1_eV = eigenvalues_eV[0]
    E2_eV = eigenvalues_eV[1]

    # Step 7: Calculate and Print the Energy Difference
    delta_E_eV = E2_eV - E1_eV

    print("The energy levels are calculated by numerically solving the Schrödinger equation.")
    print(f"The first energy level E1 is: {E1_eV:.4f} eV")
    print(f"The second energy level E2 is: {E2_eV:.4f} eV")
    print(f"\nThe energy difference is:")
    print(f"ΔE = E2 - E1 = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")
    
    return delta_E_eV

if __name__ == '__main__':
    delta_E = solve_energy_levels()
    # The final answer is enclosed in <<< >>>
    # print(f"\n<<<{delta_E:.3f}>>>")