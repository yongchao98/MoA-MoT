import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.constants import hbar, m_e, e

def solve_quantum_well_problem():
    """
    Solves the quantum well problem by numerically finding the energy levels.
    """
    # Parameters from the problem
    U0_eV = 15.0
    R_nm = 3.0
    m = m_e  # Mass of the particle (electron mass)

    # Convert parameters to SI units
    U0 = U0_eV * e
    R = R_nm * 1e-9

    # Define the corrected potential V(r) based on the plan
    def potential(r):
        """
        Calculates the potential energy V(r) in Joules based on the corrected formula.
        """
        if r < R:
            # Corrected formula for r < R: V^2(r) = U0^2 * (1 + W(exp(r/R - 1)))
            arg_exp = r / R - 1.0
            # lambertw can be complex, we take the real part of the principal branch W_0
            w_val = np.real(lambertw(np.exp(arg_exp)))
            v_sq = U0**2 * (1 + w_val)
        else:
            # Formula for r >= R: V^2(r) = U0^2 * (1 - (R/r)^2)
            # The original V0 is interpreted as U0^2 to match units.
            v_sq = U0**2 * (1 - (R / r)**2)
        
        # The potential V(r) is the square root of V^2(r).
        # The corrected formulas ensure v_sq is non-negative.
        return np.sqrt(v_sq)

    # Define the ODE system for the radial Schrödinger equation
    def ode_system(r, y, E, l):
        """
        Defines the system of first-order ODEs for the shooting method.
        y = [u, du/dr], E = energy, l = angular momentum quantum number.
        """
        u, du_dr = y
        
        # Effective potential V_eff(r) = V(r) + centrifugal term
        v_eff = potential(r)
        if r > 1e-15:  # Avoid division by zero at r=0
            v_eff += (hbar**2 * l * (l + 1)) / (2 * m * r**2)
        
        # The radial Schrödinger equation: d^2u/dr^2 = (2m/hbar^2) * (V_eff - E) * u
        d2u_dr2 = (2 * m / hbar**2) * (v_eff - E) * u
        return [du_dr, d2u_dr2]

    # Shooting method function
    def shoot(E, l, r_max_factor=10):
        """
        Solves the ODE for a given energy E and returns the value of u at r_max.
        The roots of this function are the energy eigenvalues.
        """
        r_min = 1e-15  # Start integration very close to zero
        r_max = r_max_factor * R

        # Initial conditions for u(r) near r=0, where u(r) ~ r^(l+1)
        u_min = r_min**(l + 1)
        du_dr_min = (l + 1) * r_min**l
        y0 = [u_min, du_dr_min]

        # Integrate the ODE using a standard solver
        sol = solve_ivp(
            ode_system, [r_min, r_max], y0, args=(E, l), dense_output=True
        )
        
        # Return the value of the wavefunction at the outer boundary
        return sol.sol(r_max)[0]

    # Function to find energy levels by finding roots of shoot(E, l)
    def find_energy_levels(l, num_levels, E_max_J):
        """
        Finds the first `num_levels` energy eigenvalues for a given `l` by scanning
        for roots of the shoot function.
        """
        energies = []
        E_low = 1e-6 * e  # Start search just above zero energy
        
        E1 = E_low
        psi1 = shoot(E1, l)
        
        # Scan a range of energies to find brackets where the wavefunction at r_max changes sign
        for E2 in np.linspace(E_low, E_max_J, 3000):
            psi2 = shoot(E2, l)
            if np.sign(psi1) != np.sign(psi2):
                try:
                    # A sign change indicates a root. Use brentq for robust root finding.
                    energy = brentq(shoot, E1, E2, args=(l,))
                    energies.append(energy)
                    if len(energies) == num_levels:
                        return energies
                except (ValueError, RuntimeError):
                    # Root finding might fail in some cases, just continue scanning
                    pass
            E1, psi1 = E2, psi2
        
        if not energies:
             raise RuntimeError(f"Could not find any energy levels for l={l}.")
        return energies

    # --- Main Calculation ---
    # Find the first two l=0 energy levels (n=1 and n=2 states)
    l0_energies = find_energy_levels(l=0, num_levels=2, E_max_J=U0)
    E10 = l0_energies[0]
    E20 = l0_energies[1]

    # Find the first l=1 energy level (n=1 state)
    l1_energies = find_energy_levels(l=1, num_levels=1, E_max_J=U0)
    E11 = l1_energies[0]

    # The ground state E1 is the lowest possible energy state
    E1 = E10
    # The first excited state E2 is the next lowest energy state
    E2 = min(E20, E11)

    # Calculate the final energy difference
    delta_E = E2 - E1

    # Convert results from Joules to eV for the final output
    E1_eV = E1 / e
    E2_eV = E2 / e
    delta_E_eV = delta_E / e

    # Print the results as requested
    print(f"Based on a numerical solution for the corrected potential:")
    print(f"The ground state energy is E1 = {E1_eV:.4f} eV.")
    print(f"The first excited state energy is E2 = {E2_eV:.4f} eV.")
    print(f"The energy difference is ΔE = E2 - E1 = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV.")
    
    return delta_E_eV

# Execute the solution
final_answer = solve_quantum_well_problem()
print(f"\n<<<The calculated energy difference is {final_answer:.4f} eV>>>")
