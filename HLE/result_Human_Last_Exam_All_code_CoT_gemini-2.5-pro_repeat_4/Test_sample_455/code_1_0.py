import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

def solve_quantum_well_problem():
    """
    Solves the quantum well problem to find the energy difference between the
    first two energy levels.
    """
    # --- 1. Constants and Parameters ---
    HBAR = 1.054571817e-34      # J*s
    M = 9.11e-31                # kg (mass of particle, e.g., electron)
    E_CHARGE = 1.602176634e-19  # C (electron charge for eV to Joules conversion)

    V0_eV = 15.0                # eV
    R_nm = 3.0                  # nm

    # Convert parameters to SI units
    V0_J = V0_eV * E_CHARGE
    R_m = R_nm * 1e-9

    # The potential at infinity is the limit for bound states
    V_inf_eV = np.sqrt(V0_eV)
    V_inf_J = V_inf_eV * E_CHARGE

    # Pre-calculate constant for the ODE
    ODE_CONST = 2 * M / HBAR**2

    # --- 2. Potential Function ---
    def potential_V_joules(r_m):
        """
        Calculates the potential V(r) in Joules.
        r_m: radius in meters.
        """
        # Avoid singularity at r=0 by using a small offset
        if r_m < 1e-15:
            r_m = 1e-15

        r_nm = r_m * 1e9

        if r_nm < R_nm:
            # For r < R, V^2(r) [eV^2] = V0 + W(e^(r_nm - R_nm))
            # The argument to exp must be dimensionless. We assume r and R are
            # treated as numerical values in nm.
            arg_exp = r_nm - R_nm
            # Lambert W function can return complex results, we take the real part.
            V_sq_eV = V0_eV + lambertw(np.exp(arg_exp)).real
        else:
            # For r >= R, V^2(r) [eV^2] = V0 * (1 - (R/r)^2)
            V_sq_eV = V0_eV * (1 - (R_nm / r_nm)**2)
        
        # We need sqrt(V^2), which should be positive for this problem.
        if V_sq_eV < 0:
             # This region is classically forbidden. Return a large potential.
             # This case shouldn't be hit for the bound states we seek.
            return 10 * V_inf_J
            
        V_eV = np.sqrt(V_sq_eV)
        return V_eV * E_CHARGE

    # --- 3. ODE System for Shooting Method ---
    def schrodinger_ode(r, y, E, l):
        """
        Defines the system of ODEs for the radial Schrödinger equation.
        y = [u, u'], E is energy in Joules, l is angular momentum quantum number.
        """
        u, _ = y
        
        # Effective potential V_eff = V(r) + centrifugal term
        centrifugal_term = l * (l + 1) * HBAR**2 / (2 * M * r**2)
        V_eff = potential_V_joules(r) + centrifugal_term
        
        d2u_dr2 = ODE_CONST * (V_eff - E) * u
        
        return [y[1], d2u_dr2]

    # --- 4. Function to Find Eigenvalues ---
    def find_eigenenergies(l, num_levels_to_find):
        """Finds the first `n` eigenenergies for a given l."""
        
        # Function to be used by the root-finder. It returns the value of the
        # wavefunction at a large radius r_max for a given trial energy E.
        def wavefunction_at_rmax(E, l):
            # We are looking for bound states 0 < E < V_inf
            if E <= 0 or E >= V_inf_J:
                return np.inf

            # Integration range
            r_min = 1e-15
            r_max = 20 * R_m  # Integrate to a sufficiently large radius

            # Initial conditions at r_min, based on u(r) ~ r^(l+1) for small r
            u0 = r_min**(l + 1)
            du0 = (l + 1) * r_min**l
            
            sol = solve_ivp(
                schrodinger_ode, 
                [r_min, r_max], 
                [u0, du0], 
                args=(E, l),
                dense_output=True,
                t_eval=[r_max]
            )
            
            # Return the value of u at r_max. We normalize to help the root finder.
            u_final = sol.y[0, -1]
            t_dense = np.linspace(r_min, r_max, 500)
            u_dense = sol.sol(t_dense)[0]
            max_u = np.max(np.abs(u_dense))
            return u_final / max_u if max_u > 1e-30 else u_final

        eigenenergies = []
        # Scan the allowed energy range to find brackets for the root finder
        e_scan_points = 800
        e_range = np.linspace(1e-4 * E_CHARGE, V_inf_J * 0.9999, e_scan_points)
        f_values = np.array([wavefunction_at_rmax(E, l) for E in e_range])
        
        # Find where the function crosses zero (sign change)
        sign_changes = np.where(np.diff(np.sign(f_values)))[0]
        
        for idx in sign_changes:
            if len(eigenenergies) >= num_levels_to_find:
                break
            e1, e2 = e_range[idx], e_range[idx+1]
            try:
                # Find the precise root within the bracket [e1, e2]
                root = brentq(wavefunction_at_rmax, e1, e2, args=(l,))
                eigenenergies.append(root)
            except (ValueError, RuntimeError):
                continue
        return eigenenergies

    # --- 5. Determine E1 and E2 ---
    # Find the first two s-wave (l=0) states
    energies_l0 = find_eigenenergies(l=0, num_levels_to_find=2)
    E_10_J = energies_l0[0] if len(energies_l0) > 0 else None
    E_20_J = energies_l0[1] if len(energies_l0) > 1 else None

    # Find the first p-wave (l=1) state
    energies_l1 = find_eigenenergies(l=1, num_levels_to_find=1)
    E_11_J = energies_l1[0] if len(energies_l1) > 0 else None

    if E_10_J is None or (E_20_J is None and E_11_J is None):
        print("Error: Could not find the required energy levels.")
        return

    # E1 is the ground state (lowest energy)
    E1_J = E_10_J

    # E2 is the second lowest energy state
    possible_E2 = []
    if E_20_J: possible_E2.append(E_20_J)
    if E_11_J: possible_E2.append(E_11_J)
    E2_J = min(possible_E2)

    # --- 6. Final Calculation and Output ---
    delta_E_J = E2_J - E1_J
    
    # Convert final results to eV for printing
    E1_eV = E1_J / E_CHARGE
    E2_eV = E2_J / E_CHARGE
    delta_E_eV = delta_E_J / E_CHARGE
    
    # Print the final equation with the calculated values
    print("Final Calculation:")
    print(f"First energy level (ground state) E₁ = {E1_eV:.4f} eV")
    print(f"Second energy level E₂ = {E2_eV:.4f} eV")
    print(f"ΔE = E₂ - E₁ = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")
    
    # Final answer in the required format
    print(f"\n<<<{delta_E_eV:.4f}>>>")

# Run the solver
solve_quantum_well_problem()