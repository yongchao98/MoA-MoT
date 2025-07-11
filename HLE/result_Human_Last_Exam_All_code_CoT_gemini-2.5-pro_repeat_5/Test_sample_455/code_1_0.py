import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

def solve_energy_difference():
    """
    This function calculates the energy difference between the first two energy levels
    for a particle in the specified potential well.
    """
    # --- Step 1: Define Constants and Parameters in SI units ---
    m = 9.11e-31  # Mass of the particle (kg)
    hbar = 1.054571817e-34  # Reduced Planck constant (J·s)
    e = 1.602176634e-19  # Elementary charge (C)
    V0_eV = 15.0  # Potential energy parameter (eV)
    R_nm = 3.0    # Well radius (nm)

    # Convert parameters to base SI units
    V0 = V0_eV * e
    R = R_nm * 1e-9

    # --- Step 2: Define the Corrected Potential Energy Function U(r) ---
    def potential_U(r):
        """
        Calculates the corrected potential U(r) in Joules.
        The function is vectorized to handle numpy array inputs efficiently.
        """
        @np.vectorize
        def U_scalar(r_scalar):
            if r_scalar < R:
                # This form corrects for dimensional inconsistencies in the problem statement.
                # Argument to exp is now dimensionless: r/R - 1.
                # The result is scaled by V0 to have units of energy.
                arg = np.exp(r_scalar / R - 1)
                # The principal branch (k=0) of Lambert W function is used.
                return V0 * (1 + lambertw(arg).real)
            else:  # r >= R
                if r_scalar == 0: return np.inf # Avoid division by zero
                return V0 * (1 - (R / r_scalar)**2)
        return U_scalar(r)

    def effective_potential(r, l):
        """Calculates the effective potential V_eff(r) = U(r) + centrifugal term."""
        if r == 0:
            return np.inf
        centrifugal_term = (hbar**2 * l * (l + 1)) / (2 * m * r**2)
        return potential_U(r) + centrifugal_term

    # --- Step 3: Set up the ODE and the Shooting Method ---
    def schrodinger_ode(r, y, E, l):
        """Defines the ODE system for the radial Schrödinger equation."""
        u, du_dr = y
        V_eff = effective_potential(r, l)
        d2u_dr2 = (2 * m / hbar**2) * (V_eff - E) * u
        return [du_dr, d2u_dr2]

    # Integration range for the ODE solver
    r_min = 1e-15  # Start integration slightly away from r=0 to avoid singularity
    r_max = 10 * R  # Integrate far enough for the bound state wavefunction to decay

    def get_wavefunction_end_value(E, l):
        """
        Solves the ODE for a given energy E and returns the value of u(r_max).
        We will find the roots of this function to get the energy eigenvalues.
        """
        y0 = [0, 1]  # Initial conditions: u(r_min)=0, u'(r_min)=1
        sol = solve_ivp(
            fun=schrodinger_ode, t_span=[r_min, r_max], y0=y0, args=(E, l), dense_output=True
        )
        u_at_rmax = sol.sol(r_max)[0]
        return u_at_rmax

    # --- Step 4: Find the Energy Eigenvalues E1 and E2 ---
    def find_nth_energy_level(n, l):
        """
        Finds the n-th energy level for a given angular momentum l by finding
        brackets [a,b] where the shooting function changes sign, and then using
        a root-finder (brentq) within those brackets.
        """
        # Lower bound for energy search is the minimum of the effective potential.
        # Bound states must have energy E < V0.
        res = minimize_scalar(lambda r: effective_potential(r,l), bounds=(r_min, r_max), method='bounded')
        E_min_search = res.fun
        E_max_search = V0

        found_energies = []
        scan_points = np.linspace(E_min_search, E_max_search, 2000)
        f_prev = get_wavefunction_end_value(scan_points[0], l)

        for i in range(1, len(scan_points)):
            E_prev, E_curr = scan_points[i-1], scan_points[i]
            f_curr = get_wavefunction_end_value(E_curr, l)
            
            # If sign changes, a root is in the interval [E_prev, E_curr]
            if np.sign(f_prev) != np.sign(f_curr):
                try:
                    root = brentq(lambda E: get_wavefunction_end_value(E, l), E_prev, E_curr)
                    found_energies.append(root)
                    if len(found_energies) == n:
                        return found_energies[n-1]
                except (RuntimeError, ValueError):
                    # brentq might fail in some edge cases, continue scanning
                    pass
            f_prev = f_curr

        if len(found_energies) < n:
            raise RuntimeError(f"Could not find {n} energy levels for l={l}.")
        return found_energies[n-1]

    # E1 is the ground state (n=1, l=0)
    E1_J = find_nth_energy_level(n=1, l=0)

    # E2 is the first excited state, which is min(E_n=2,l=0, E_n=1,l=1)
    E_2_0_J = find_nth_energy_level(n=2, l=0)
    E_1_1_J = find_nth_energy_level(n=1, l=1)
    E2_J = min(E_2_0_J, E_1_1_J)

    # --- Step 5: Calculate the Energy Difference and Print ---
    delta_E_J = E2_J - E1_J

    # Convert final values from Joules to electron-Volts (eV) for the output
    E1_eV = E1_J / e
    E2_eV = E2_J / e
    delta_E_eV = delta_E_J / e

    print(f"The first energy level (ground state) is E1 = {E1_eV:.4f} eV.")
    print(f"The second energy level (first excited state) is E2 = {E2_eV:.4f} eV.")
    print("\nThe energy difference is:")
    print(f"ΔE = E2 - E1 = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")
    
    # Return the final numerical answer for the submission format
    return delta_E_eV

# Execute the calculation and print the results
final_answer = solve_energy_difference()
print(f"\n<<<final_answer>>>\n{final_answer:.4f}")