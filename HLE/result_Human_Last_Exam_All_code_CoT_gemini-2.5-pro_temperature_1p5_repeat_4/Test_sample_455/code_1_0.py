import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_energy_levels():
    """
    Calculates the energy difference between the first and second energy levels
    of a quantum particle in a given potential well.
    """
    # --- 1. Define Constants and Potential ---

    # Parameters from the problem
    V0_sq = 15.0  # V_0 in the formula, units are effectively eV^2
    R = 3.0       # Radius in nm
    # Mass of the particle (electron) m = 9.11e-31 kg
    
    # Physical constant hbar^2 / (2 * m_e) in units of eV * nm^2
    HBAR2_2M = 0.0380998

    def potential_V(r):
        """
        Calculates the potential V(r) in eV, given r in nm.
        The input formula is for V^2(r).
        """
        if r < R:
            # For r < R, V^2(r) = V0 + W(exp(r - R))
            # lambertw requires a complex argument and can return a complex result.
            # We take the real part as the physical potential is real.
            # For exp(r-R) with r<R, the argument to W is real and < 1.
            val_V2 = V0_sq + lambertw(np.exp(r - R)).real
        else:
            # For r >= R, V^2(r) = V0 * (1 - (r/R)^-2)
            if r == 0:  # Should not happen in this branch, but for safety
                return np.inf
            # Clip is used to prevent taking a sqrt of a tiny negative number 
            # due to floating point inaccuracies near r = R.
            term = np.clip(1.0 - (R / r)**2, 0, None)
            val_V2 = V0_sq * term
        
        # Potential should not be negative. V^2 must be >= 0.
        if val_V2 < 0:
            return np.inf
        return np.sqrt(val_V2)

    def effective_potential_Veff(r, l):
        """Calculates the effective potential including the centrifugal term."""
        if r == 0:
            return np.inf
        return potential_V(r) + HBAR2_2M * l * (l + 1) / r**2

    # --- 2. Implement the Shooting Method ---

    def schrodinger_ode(r, y, E, l):
        """
        Defines the system of first-order ODEs for the radial Schrodinger equation.
        y = [u, u']
        returns dy/dr = [u', u'']
        """
        u, _ = y
        # From -C*u'' + V_eff*u = E*u  =>  u'' = (V_eff - E)*u / C
        d2u_dr2 = (effective_potential_Veff(r, l) - E) * u / HBAR2_2M
        return [y[1], d2u_dr2]

    def shoot(E, l, r_max=20.0, r_min=1e-6):
        """
        Integrates the ODE for a given energy E and returns the wavefunction
        value at the outer boundary r_max. The roots of this function are the eigenvalues.
        """
        # Start integration close to zero with initial conditions u(0)=0 and u'(0)=1
        y0 = [0.0, 1.0]
        sol = solve_ivp(
            schrodinger_ode,
            [r_min, r_max],
            y0,
            args=(E, l)
        )
        # Return the value of u(r) at the last integration point
        return sol.y[0, -1]

    def find_eigenvalues(l, num_levels, bracket_range):
        """Finds the first 'num_levels' eigenvalues for a given 'l'."""
        E_min, E_max = bracket_range
        energies = []
        
        # Scan for brackets where the shooting function crosses zero
        search_points = np.linspace(E_min, E_max, 300)
        shoot_values = np.array([shoot(E, l) for E in search_points])
        
        # Find indices where sign changes occur
        sign_changes = np.where(np.sign(shoot_values[:-1]) != np.sign(shoot_values[1:]))[0]
        
        for i in sign_changes:
            if len(energies) >= num_levels:
                break
            e1, e2 = search_points[i], search_points[i+1]
            try:
                # Find the root within the bracket
                eigen_E = brentq(shoot, e1, e2, args=(l,))
                energies.append(eigen_E)
            except ValueError:
                continue
        return energies

    # --- 3. Calculate Energy Levels and the Difference ---

    # The potential at infinity is sqrt(V0) ~= 3.87 eV. Bound states must be below this.
    search_range = (0.01, np.sqrt(V0_sq) - 0.01)
    
    # Find the two lowest s-states (l=0)
    e_l0 = find_eigenvalues(l=0, num_levels=2, bracket_range=search_range)
    E_10 = e_l0[0]  # Ground state (n=1, l=0)
    E_20 = e_l0[1]  # First excited s-state (n=2, l=0)

    # Find the lowest p-state (l=1)
    # The lowest l=1 state must have higher energy than the l=0 ground state.
    search_range_l1 = (E_10, search_range[1])
    e_l1 = find_eigenvalues(l=1, num_levels=1, bracket_range=search_range_l1)
    E_11 = e_l1[0]  # Ground p-state (n=1, l=1)
    
    # The first energy level (ground state) is E1
    E1 = E_10
    
    # The second energy level is the minimum of the next available states
    E2 = min(E_20, E_11)
    
    # The energy difference
    Delta_E = E2 - E1

    # --- 4. Print the Final Result ---
    
    print(f"The ground state energy is E1 = {E1:.4f} eV.")
    print(f"The second energy level is E2 = {E2:.4f} eV.")
    print(f"The energy difference is:")
    print(f"ΔE = E2 - E1 = {E2:.4f} eV - {E1:.4f} eV = {Delta_E:.4f} eV")
    
    return Delta_E

# Execute the calculation and capture the final numeric answer
final_answer = solve_energy_levels()

# The final answer in the required format
print(f"\n<<<ΔE = {final_answer:.4f} eV>>>")