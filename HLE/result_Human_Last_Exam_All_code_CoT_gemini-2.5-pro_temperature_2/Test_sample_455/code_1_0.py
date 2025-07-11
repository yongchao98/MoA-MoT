import numpy as np
from scipy.special import lambertw
from scipy.integrate import odeint
from scipy.optimize import brentq
import scipy.constants as const

def solve_quantum_well_energy_difference():
    """
    Calculates the energy difference between the first two energy levels of a particle
    in the specified potential well.
    """
    # --- Setup and Constants ---
    V0_eV_val = 15.0  # Numerical value for V0 from the problem
    R_nm = 3.0
    R = R_nm * 1e-9  # Convert radius to meters
    m = const.m_e
    hbar = const.hbar
    e = const.e
    const_factor_schrodinger = 2 * m / hbar**2

    # --- Potential Function ---
    # Interprets the given formula for V^2(r) as yielding a value in eV^2.
    def potential_in_joules(r):
        """Calculates potential V(r) in Joules for a given radius r in meters."""
        # Use a small floor for r to avoid division by zero or other numerical issues.
        r = max(1e-15, r)
        r_nm = r * 1e9

        if r_nm < R_nm:
            arg_w = np.exp(r_nm - R_nm)
            # Use the real part of the principal branch of the Lambert W function.
            w_val = np.real(lambertw(arg_w))
            v_sq_ev = V0_eV_val + w_val
        else:
            v_sq_ev = V0_eV_val * (1.0 - (R_nm / r_nm)**2)

        # The square of the potential should be non-negative for a real potential.
        v_ev = np.sqrt(max(0, v_sq_ev))
        return v_ev * e

    # --- Schrödinger Equation ODE System for the Solver ---
    def radial_ode(y, r, E_J, l):
        """Defines the ODE system for u(r) based on the radial Schrödinger equation."""
        u, du_dr = y
        V_eff = potential_in_joules(r) + (hbar**2 * l * (l + 1)) / (2 * m * r**2)
        d2u_dr2 = const_factor_schrodinger * (V_eff - E_J) * u
        return [du_dr, d2u_dr2]

    # --- Shooting Method Function ---
    def find_wavefunction_at_rmax(E_J, l, r_max, num_points):
        """Solves the radial ODE for energy E and returns u(r_max)."""
        r_min = 1e-15
        r_span = np.linspace(r_min, r_max, num_points)

        # Set initial conditions for u(r) near r=0 based on l.
        if l == 0:
            u0 = [r_min, 1.0]  # For s-wave, u(r) is proportional to r
        else:
            u0 = [r_min**(l + 1), (l + 1) * r_min**l] # u(r) ~ r^(l+1)

        solution = odeint(radial_ode, u0, r_span, args=(E_J, l))
        return solution[-1, 0] # Return the value of u at r_max

    # --- Function to Find Energy Eigenvalues ---
    def find_energy_levels(l, n_levels, E_min_eV, E_max_eV):
        """Scans an energy range to find the first 'n_levels' for a given l."""
        r_max = 8 * R  # Integrate to a sufficiently large radius
        num_points = 2000
        energies_eV = []
        
        # Scan a range of energies to find brackets where the wavefunction solution crosses zero.
        E_scan_eV = np.linspace(E_min_eV, E_max_eV, 300)
        u_values = [find_wavefunction_at_rmax(E * e, l, r_max, num_points) for E in E_scan_eV]

        # Find where sign changes occur, indicating a root.
        sign_changes = np.where(np.diff(np.sign(u_values)))[0]

        for i in sign_changes:
            if len(energies_eV) >= n_levels:
                break
            # Use a precise root-finder (Brent's method) within the found bracket.
            e_low_eV, e_high_eV = E_scan_eV[i], E_scan_eV[i+1]
            try:
                eigen_energy_J = brentq(find_wavefunction_at_rmax, e_low_eV * e, e_high_eV * e, args=(l, r_max, num_points))
                energies_eV.append(eigen_energy_J / e)
            except ValueError:
                # This can happen if the sign at the boundaries is not different; ignore.
                continue
                
        return energies_eV

    # The potential's minimum is 0 eV, and it approaches sqrt(15) ~= 3.87 eV at infinity.
    # Bound states must exist within this energy range.
    E_min_search_eV = 0.1
    E_max_search_eV = 3.85

    # Find E1: The ground state energy is the lowest l=0 state (1s state).
    E_levels_l0 = find_energy_levels(l=0, n_levels=2, E_min_eV=E_min_search_eV, E_max_eV=E_max_search_eV)
    if not E_levels_l0:
        print("Could not find the ground state energy level (l=0).")
        return
    E1 = E_levels_l0[0]

    # Find E2: The first excited state is the minimum of the second l=0 state (2s) 
    # and the first l=1 state (1p).
    if len(E_levels_l0) < 2:
        print("Could not find the second l=0 energy level.")
        return
    E2_candidate1_2s = E_levels_l0[1]

    E_levels_l1 = find_energy_levels(l=1, n_levels=1, E_min_eV=E_min_search_eV, E_max_eV=E_max_search_eV)
    if not E_levels_l1:
        print("Could not find the first l=1 energy level.")
        return
    E2_candidate2_1p = E_levels_l1[0]
    
    E2 = min(E2_candidate1_2s, E2_candidate2_1p)
    
    # Calculate the energy difference.
    delta_E_eV = E2 - E1
    
    # Output the final result as an equation.
    print("The first energy level E1 (1s state) was found at:")
    print(f"E1 = {E1:.4f} eV")
    print("\nThe second energy level E2 is the minimum of the 2s and 1p states:")
    print(f"E(2s) = {E2_candidate1_2s:.4f} eV")
    print(f"E(1p) = {E2_candidate2_1p:.4f} eV")
    print(f"E2 = min(E(2s), E(1p)) = {E2:.4f} eV")
    print("\nThe energy difference is Delta E = E2 - E1:")
    print(f"{E2:.4f} - {E1:.4f} = {delta_E_eV:.4f}")

solve_quantum_well_energy_difference()
<<<0.9080>>>