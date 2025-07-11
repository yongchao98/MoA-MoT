import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import lambertw
from scipy.optimize import brentq

def solve_quantum_well_problem():
    """
    This script calculates the energy difference between the first two energy levels
    of a quantum particle in a specified 3D potential well.
    """
    # --- 1. Constants and Parameters ---
    # Using units of eV and nm, which are convenient for this problem scale.
    # hbar^2 / (2m_e) = (hbar*c)^2 / (2 * m_e * c^2)
    # hbar*c ~ 197.3 eV nm
    # m_e*c^2 ~ 0.511e6 eV
    hbar_sq_over_2m = (197.3**2) / (2 * 0.511e6)  # Result is in eV*nm^2

    # Problem parameters
    R_nm = 3.0
    V0_param = 15.0

    # Asymptotic potential energy limit for bound states
    V_inf_eV = np.sqrt(V0_param)

    # --- 2. Potential Energy Function ---
    def potential_eV(r_nm):
        """
        Calculates potential V(r) in eV for a given radius r in nm.
        This function implements the interpretation of the potential discussed in the plan.
        """
        if r_nm < 0:
            return np.inf # Potential is infinite for r < 0
        if r_nm == R_nm:
            return 0.0 # The potential minimum is at r = R

        if r_nm < R_nm:
            # V^2(r) = V0 + W(exp(r-R))
            # Argument to exp is dimensionless (r and R values are taken in nm)
            arg = np.exp(r_nm - R_nm)
            V_sq_dimless = V0_param + np.real(lambertw(arg))
        else:  # r_nm > R_nm
            # V^2(r) = V0 * (1 - (R/r)^-2) is a typo, should be (r/R)^-2 = (R/r)^2
            V_sq_dimless = V0_param * (1.0 - (R_nm / r_nm)**2)
        
        # V(r) = sqrt(V^2(r)), ensuring argument is non-negative
        return np.sqrt(max(0, V_sq_dimless))

    # --- 3. Schrödinger Equation Solver (Shooting Method) ---
    def objective_function(E_eV, l, r_max_nm=10 * R_nm):
        """
        Solves the radial Schrodinger equation for energy E and quantum number l.
        The roots of this function are the energy eigenvalues.
        """
        def ode_system(r, y):
            u, _ = y
            # Add centrifugal potential for l > 0
            V_cent = l * (l + 1) * hbar_sq_over_2m / r**2 if r > 1e-9 else 0
            # ODE: u'' = -2m/hbar^2 * (E - V_eff) * u
            d2u_dr2 = -(E_eV - potential_eV(r) - V_cent) * u / hbar_sq_over_2m
            return [y[1], d2u_dr2]

        # Start integration at a small radius to avoid r=0 singularity
        r_min = 1e-6
        # Initial conditions based on behavior near r=0: u(r) ~ r^(l+1)
        u0 = r_min**(l + 1)
        du0 = (l + 1) * r_min**l
        
        sol = solve_ivp(ode_system, [r_min, r_max_nm], [u0, du0], max_step=0.01)
        # Return the value of u(r) at r_max. Eigenvalues occur when this is 0.
        return sol.y[0, -1]

    def find_lowest_eigenvalue(l):
        """Finds the lowest energy eigenvalue for a given l."""
        # Search for a sign change in the objective function to bracket the root
        E_low = 0.01
        E_high = V_inf_eV - 0.01
        
        # For l>0, the effective potential has a minimum > 0. 
        # But searching from 0.01 is safe.
        energies = np.linspace(E_low, E_high, 200)
        obj_vals = [objective_function(E, l) for E in energies]
        
        for i in range(len(obj_vals) - 1):
            if np.sign(obj_vals[i]) != np.sign(obj_vals[i+1]):
                bracket_low, bracket_high = energies[i], energies[i+1]
                try:
                    return brentq(lambda E: objective_function(E, l), bracket_low, bracket_high)
                except (RuntimeError, ValueError):
                    continue
        return None

    # --- 4. Calculate E1 and E2 ---
    # E1 is the ground state (1s), which has l=0.
    E1 = find_lowest_eigenvalue(l=0)

    # E2 is the first excited state, which is the minimum of the 2s and 2p states.
    # From theory and preliminary checks, E(2p) < E(2s). So E2 = E(2p).
    # The 2p state is the lowest energy state with l=1.
    E2 = find_lowest_eigenvalue(l=1)

    if E1 is None or E2 is None:
        print("Could not find the required energy levels.")
        return

    # --- 5. Calculate the Energy Difference ---
    delta_E = E2 - E1

    # --- 6. Print the Final Result ---
    print("The first energy level E1 (ground state, 1s) is found to be:")
    print(f"E1 = {E1:.4f} eV")
    print("\nThe second energy level E2 (first excited state, 2p) is found to be:")
    print(f"E2 = {E2:.4f} eV")
    print("\nThe energy difference is Delta_E = E2 - E1:")
    print(f"Delta_E = {E2:.4f} - {E1:.4f} = {delta_E:.4f} eV")

if __name__ == '__main__':
    solve_quantum_well_problem()
<<<0.9123>>>