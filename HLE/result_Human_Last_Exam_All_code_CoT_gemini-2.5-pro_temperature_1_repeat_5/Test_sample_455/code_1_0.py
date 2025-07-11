import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_quantum_well():
    """
    Solves the Schrodinger equation for the given potential to find the
    energy difference between the first two energy levels.
    """
    # --- Step 1: Define constants and parameters ---
    # Physical constants in SI units
    m_particle = 9.11e-31  # kg (mass of electron, given as particle mass)
    hbar = 1.054571817e-34 # J*s
    e_charge = 1.602176634e-19 # C

    # Problem parameters based on interpretation
    V0_param_eV2 = 15.0  # Interpreted as 15 eV^2
    R_param_nm = 3.0     # Radius in nm

    # Convert parameters to SI units for calculations
    R_m = R_param_nm * 1e-9  # meters

    # Pre-calculate constant for ODE: 2*m/hbar^2
    C1 = 2 * m_particle / (hbar**2) # 1 / (J * m^2)

    print("This script solves for the energy levels in the given potential well.")
    print("Assumptions made to resolve dimensional inconsistencies in the potential function:")
    print(f"1. V_0 = {V0_param_eV2} is interpreted as being in units of eV^2.")
    print(f"2. In the term W(exp(r-R)), r and R are the numerical values in nanometers.")
    print("-" * 20)

    # --- Step 2: Define the potential energy function ---
    def potential_squared_J2(r_m):
        """
        Calculates the square of the potential energy (V^2) in Joules^2.
        Input r is in meters.
        """
        r_nm = r_m * 1e9 # Convert radius to nanometers for the formula

        if np.isclose(r_nm, R_param_nm):
            # The potential is discontinuous at r=R. This is the well bottom (limit from r > R).
            return 0.0

        if 0 <= r_nm < R_param_nm:
            # For r < R: V^2 [eV^2] = 15 + W(e^(r_nm-R_nm))
            V_sq_eV2 = V0_param_eV2 + np.real(lambertw(np.exp(r_nm - R_param_nm)))
            return V_sq_eV2 * (e_charge**2)
        else: # r_nm > R_param_nm
            # For r >= R: V^2 [eV^2] = 15 * (1 - (r_nm/R_nm)^-2)
            V_sq_eV2 = V0_param_eV2 * (1 - (r_nm / R_param_nm)**(-2))
            return V_sq_eV2 * (e_charge**2)

    def effective_potential_J(r_m, l):
        """
        Calculates the effective potential in Joules, including the centrifugal term.
        """
        if r_m == 0:
            return np.inf if l > 0 else np.sqrt(potential_squared_J2(0))

        V_r = np.sqrt(potential_squared_J2(r_m))
        centrifugal_term = l * (l + 1) * (hbar**2) / (2 * m_particle * r_m**2)
        return V_r + centrifugal_term

    # --- Step 3: Set up the shooting method ---
    r_min = 1e-14 # Start integration close to zero (in meters)
    r_max = 20 * R_m # Integrate out to a large radius

    def schrodinger_ode(r, y, energy_J, l):
        """Defines the ODE system: y' = [u', u'']"""
        u, _ = y
        V_eff_J = effective_potential_J(r, l)
        d2u_dr2 = -C1 * (energy_J - V_eff_J) * u
        return [y[1], d2u_dr2]

    def shoot(energy_J, l):
        """
        Solves the ODE for a given energy E and returns u(r_max).
        Eigenvalues are energies E for which u(r_max) -> 0.
        """
        # Boundary conditions near r=0: u(r) ~ r^(l+1)
        y0 = [0.0, 1.0] if l == 0 else [r_min**(l+1), (l+1) * r_min**l]

        sol = solve_ivp(
            fun=schrodinger_ode, t_span=[r_min, r_max], y0=y0,
            args=(energy_J, l), method='RK45'
        )
        return sol.y[0, -1]

    # --- Step 4: Find the energy levels by finding roots ---
    def find_energy_levels(l, num_levels_to_find):
        """Finds the first `num_levels` eigenvalues for a given l."""
        levels_J = []
        V_inf_J = np.sqrt(V0_param_eV2 * e_charge**2)
        E_min_J = 1e-22 # Start search just above zero
        E_max_J = V_inf_J * 0.999 # Search up to the potential at infinity

        # Scan the energy range to bracket the roots
        e_scan = np.linspace(E_min_J, E_max_J, 750)
        u_final = np.array([shoot(E, l) for E in e_scan])
        
        sign_changes = np.where(np.diff(np.sign(u_final)))[0]
        
        for idx in sign_changes:
            if len(levels_J) >= num_levels_to_find: break
            E1_bracket, E2_bracket = e_scan[idx], e_scan[idx+1]
            try:
                eigenvalue_J = brentq(shoot, E1_bracket, E2_bracket, args=(l,))
                levels_J.append(eigenvalue_J)
            except (ValueError, RuntimeError):
                continue
        return np.array(levels_J)

    # --- Step 5: Calculate the final result ---
    print("Searching for the lowest energy states (l=0 and l=1)...")
    # Find the first two l=0 energy levels (n=1, n=2)
    levels_l0_J = find_energy_levels(l=0, num_levels_to_find=2)
    E10_J = levels_l0_J[0] if len(levels_l0_J) > 0 else np.inf
    E20_J = levels_l0_J[1] if len(levels_l0_J) > 1 else np.inf

    # Find the first l=1 energy level (n=1)
    levels_l1_J = find_energy_levels(l=1, num_levels_to_find=1)
    E11_J = levels_l1_J[0] if len(levels_l1_J) > 0 else np.inf

    all_levels_J = sorted([E10_J, E20_J, E11_J])

    E1_J = all_levels_J[0]
    E2_J = all_levels_J[1]

    # Convert final results to eV
    E1_eV = E1_J / e_charge
    E2_eV = E2_J / e_charge
    delta_E_eV = (E2_J - E1_J) / e_charge
    
    print("-" * 20)
    print(f"The ground state energy (E1) is found to be {E1_eV:.4f} eV.")
    print(f"The first excited state energy (E2) is found to be {E2_eV:.4f} eV.")
    print("\nThe energy difference is Delta E = E2 - E1.")
    print(f"Final Equation: {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")
    
    return delta_E_eV

if __name__ == '__main__':
    final_answer = solve_quantum_well()
    # The required format for the final answer
    # print(f"<<<{final_answer:.4f}>>>")