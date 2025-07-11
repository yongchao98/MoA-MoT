import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import lambertw
from scipy.constants import hbar, m_e, e

def solve_quantum_well_problem():
    """
    This function encapsulates the full solution to the quantum well problem.
    It calculates the energy difference between the first and second energy levels.
    """
    # --- Parameters and Constants ---
    m = m_e             # Mass of the particle (electron mass)
    R_nm = 3.0          # Radius of the well in nm
    R = R_nm * 1e-9     # Radius in SI units (meters)
    V0_param = 15.0     # The 'V0' parameter from the formula, value is 15.

    # --- Potential Energy Function ---
    def potential_in_eV(r_m):
        """
        Calculates potential V(r) in eV for a radius r in meters.
        """
        r_arr = np.atleast_1d(r_m)
        V_sq_eV2 = np.zeros_like(r_arr, dtype=float)

        # For 0 <= r < R
        mask1 = (r_arr >= 0) & (r_arr < R)
        r1 = r_arr[mask1]
        if r1.size > 0:
            arg = (r1 - R) / R
            # The lambertw function can return complex results for negative args,
            # but our argument exp(...) is always positive. We take the real part for safety.
            lambertw_val = np.real(lambertw(np.exp(arg)))
            V_sq_eV2[mask1] = V0_param + lambertw_val

        # For r >= R
        mask2 = r_arr >= R
        r2 = r_arr[mask2]
        if r2.size > 0:
            V_sq_eV2[mask2] = V0_param * (1 - (r2 / R)**(-2))

        potential_eV = np.sqrt(V_sq_eV2)
        return potential_eV[0] if isinstance(r_m, (int, float)) else potential_eV

    # --- Schrödinger ODE System ---
    def schrodinger_ode_sys(r, y, E_J, l):
        """
        System of ODEs for the radial Schrödinger equation.
        y=[u, u'], E_J is energy in Joules, l is angular momentum number.
        """
        u, du_dr = y
        # The ODE solver does not evaluate exactly at r=0
        if r == 0: return [0, 0]

        V_J = potential_in_eV(r) * e
        centrifugal_J = (hbar**2 * l * (l + 1)) / (2 * m * r**2) if l > 0 else 0
        V_eff_J = V_J + centrifugal_J
        d2u_dr2 = (2 * m / hbar**2) * (V_eff_J - E_J) * u
        return [du_dr, d2u_dr2]

    # --- Eigenenergy Solver using the Shooting Method ---
    def find_eigen_energy_joules(l, num_radial_nodes):
        """
        Finds the eigenenergy in Joules for a state (l, num_radial_nodes).
        """
        r_max = 12 * R      # Integration limit, far enough for wavefunction to decay
        r_small = 1e-15     # Starting point to avoid r=0 singularity

        # Initial conditions for u(r) and u'(r) near r=0
        y0 = [r_small, 1.0] if l == 0 else [r_small**(l + 1), (l + 1) * r_small**l]

        def objective_function(E_J):
            """Returns the wavefunction value at r_max for a given energy."""
            sol = solve_ivp(schrodinger_ode_sys, [r_small, r_max], y0, args=(E_J, l), method='RK45')
            return sol.y[0, -1]

        # Search for energy roots in a plausible range (0 to the potential barrier height at infinity)
        E_max_eV = np.sqrt(V0_param)
        E_min_J, E_max_J = 0.01 * e, (E_max_eV - 0.01) * e

        energies_to_scan = np.linspace(E_min_J, E_max_J, 250)
        obj_values = [objective_function(E) for E in energies_to_scan]

        # Find energy brackets where the wavefunction at r_max changes sign
        eigenenergies = []
        for i in range(len(obj_values) - 1):
            if np.sign(obj_values[i]) != np.sign(obj_values[i+1]):
                try:
                    # Find the precise root (eigenenergy) within the bracket
                    root = brentq(objective_function, energies_to_scan[i], energies_to_scan[i+1])
                    eigenenergies.append(root)
                except (ValueError, RuntimeError):
                    continue

        # For each found eigenenergy, count its nodes to identify the correct state
        for E_J in sorted(eigenenergies):
            sol = solve_ivp(schrodinger_ode_sys, [r_small, r_max], y0, args=(E_J, l), dense_output=True)
            r_eval = np.linspace(r_small, r_max, 2000)
            u_r = sol.sol(r_eval)[0]
            nodes = len(np.where(np.diff(np.sign(u_r)))[0])
            if nodes == num_radial_nodes:
                return E_J
        return None

    # --- Main Calculation ---
    # E1: Ground state (1s -> l=0, 0 radial nodes)
    E1_J = find_eigen_energy_joules(l=0, num_radial_nodes=0)
    E1_eV = E1_J / e

    # E(2s): l=0, 1 radial node
    E_2s_J = find_eigen_energy_joules(l=0, num_radial_nodes=1)
    E_2s_eV = E_2s_J / e

    # E(1p): l=1, 0 radial nodes
    E_1p_J = find_eigen_energy_joules(l=1, num_radial_nodes=0)
    E_1p_eV = E_1p_J / e

    # E2: First excited state energy is the minimum of the two possibilities
    E2_eV = min(E_2s_eV, E_1p_eV)

    # Final energy difference
    delta_E_eV = E2_eV - E1_eV

    # --- Print Final Answer ---
    print(f"Ground state energy (E1): {E1_eV:.4f} eV")
    print(f"First excited state energy (E2): {E2_eV:.4f} eV")
    print("\nFinal Calculation in eV:")
    print(f"ΔE = E2 - E1 = {E2_eV:.4f} - {E1_eV:.4f} = {delta_E_eV:.4f}")

    # Output the final numerical answer in the required format
    global final_answer
    final_answer = f"<<<{delta_E_eV:.4f}>>>"


if __name__ == '__main__':
    final_answer = ""
    solve_quantum_well_problem()
    print(f"\n{final_answer}")