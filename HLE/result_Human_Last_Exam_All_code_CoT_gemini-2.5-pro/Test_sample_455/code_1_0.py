import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_quantum_well_energy_difference():
    """
    Calculates the energy difference between the first two energy levels of a particle
    in a given 3D potential well using a numerical shooting method.
    """
    # --- 1. Constants and Parameters (in SI units) ---
    V0_eV = 15.0
    R_nm = 3.0
    m = 9.11e-31  # mass of the particle (electron) in kg
    e_charge = 1.602e-19 # elementary charge in C
    hbar = 1.05457e-34 # reduced Planck constant in J*s

    # Convert parameters to SI units
    V0 = V0_eV * e_charge
    R = R_nm * 1e-9

    # --- 2. Simplified Potential and Schrödinger Equation ---
    # Based on the analysis, we solve for r >= R with a hard wall at r=R.
    def potential_outer(r):
        """Potential V(r) for r >= R. Returns potential in Joules."""
        if r < R:
            return np.inf # This region is inaccessible
        arg = 1.0 - (R / r)**2
        return np.sqrt(V0 * arg)

    def schrodinger_ode(r, y, E):
        """
        Defines the system of first-order ODEs for the radial Schrödinger equation:
        y = [u, u'], dy/dr = [u', u'']
        """
        u, du_dr = y
        V_r = potential_outer(r)
        d2u_dr2 = (2 * m / hbar**2) * (V_r - E) * u
        return [du_dr, d2u_dr2]

    # --- 3. Backward Shooting Method Implementation ---
    def find_wavefunction_at_R(E):
        """
        Solves the ODE for a given energy E by integrating backwards from a large r
        down to r=R. It returns the value of the wavefunction u(R).
        The energy eigenvalues are the roots of this function.
        """
        V_inf = np.sqrt(V0)
        # A bound state must have energy 0 < E < V_inf
        if not (0 < E < V_inf):
            return 1e10 # Return a large number if not a valid energy

        # Start integration far enough for the asymptotic solution to be valid
        r_start = 20 * R
        r_end = R
        
        # Asymptotic behavior for u(r) at large r is u(r) ~ exp(-kappa*r)
        # where kappa = sqrt(2*m*(V_inf - E))/hbar
        kappa = np.sqrt(2 * m * (V_inf - E)) / hbar
        
        # Set initial conditions at r_start. Normalization is arbitrary.
        u_start = np.exp(-kappa * r_start)
        du_dr_start = -kappa * u_start
        initial_conditions = [u_start, du_dr_start]
        
        # Solve the ODE Initial Value Problem
        sol = solve_ivp(
            fun=schrodinger_ode,
            t_span=[r_start, r_end], # Integrate backwards in r
            y0=initial_conditions,
            args=(E,),
            method='DOP853', # A precise, adaptive solver
            rtol=1e-9,
            atol=1e-9
        )
        
        # Return the wavefunction's value at the endpoint r=R
        return sol.y[0, -1]

    # --- 4. Find Energy Eigenvalues ---
    eigenvalues_J = []
    V_inf_eV = np.sqrt(V0_eV)
    
    # Scan for energies where u(R) crosses zero, indicating an eigenvalue
    energy_scan_eV = np.linspace(0.01, V_inf_eV - 0.01, 1000)
    u_at_R_values = [find_wavefunction_at_R(E * e_charge) for E in energy_scan_eV]
    
    # Find the roots (eigenvalues) using the brentq algorithm for precision
    for i in range(len(energy_scan_eV) - 1):
        if np.sign(u_at_R_values[i]) != np.sign(u_at_R_values[i+1]):
            E_low_J = energy_scan_eV[i] * e_charge
            E_high_J = energy_scan_eV[i+1] * e_charge
            
            try:
                eigenvalue_J = brentq(find_wavefunction_at_R, E_low_J, E_high_J)
                eigenvalues_J.append(eigenvalue_J)
                if len(eigenvalues_J) >= 2: # Stop after finding the first two
                    break
            except ValueError:
                # brentq fails if the signs are the same, which shouldn't happen here
                pass
    
    if len(eigenvalues_J) < 2:
        print("Error: Could not find two energy levels in the specified range.")
        return

    # --- 5. Calculate and Print the Final Result ---
    E1_J, E2_J = eigenvalues_J[0], eigenvalues_J[1]

    # Convert final energies from Joules to eV
    E1_eV = E1_J / e_charge
    E2_eV = E2_J / e_charge

    # Calculate the energy difference
    delta_E_eV = E2_eV - E1_eV
    
    print(f"The first energy level E1 is {E1_eV:.4f} eV.")
    print(f"The second energy level E2 is {E2_eV:.4f} eV.")
    print(f"The energy difference is ΔE = E2 - E1 = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV.")

solve_quantum_well_energy_difference()