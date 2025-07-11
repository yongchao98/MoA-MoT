import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import lambertw
from scipy.optimize import brentq

def solve_quantum_well_energy_difference():
    """
    Calculates the energy difference between the first and second energy levels
    for a particle in the given 3D potential well.
    """
    # 1. Define physical constants and problem parameters
    
    # Mass of the particle (electron) in kg
    m = 9.11e-31
    # Planck's constant (h-bar) in J.s
    hbar = 1.054571817e-34
    # Conversion factor from eV to Joules
    eV_to_J = 1.602176634e-19

    # Pre-calculate hbar^2 / (2m) in units of eV * nm^2
    hbar_sq_2m_SI = hbar**2 / (2 * m)
    hbar_sq_2m = hbar_sq_2m_SI / eV_to_J * 1e18

    # Potential parameters
    V0 = 15.0  # eV
    R = 3.0    # nm

    # 2. Implement the potential energy function V(r)
    
    def potential(r):
        """
        Calculates the potential V(r) in eV at a given radius r in nm.
        This function is vectorized to handle numpy arrays efficiently.
        """
        # Ensure r is a numpy array for vectorized operations
        if np.isscalar(r):
            is_scalar = True
            r = np.array([r])
        else:
            is_scalar = False
            
        v_sq = np.zeros_like(r, dtype=float)
        
        # Region 1: 0 <= r < R
        mask_lt_R = r < R
        if np.any(mask_lt_R):
            r_lt_R = r[mask_lt_R]
            v_sq[mask_lt_R] = V0 + lambertw(np.exp(r_lt_R - R)).real
        
        # Region 2: r >= R
        mask_ge_R = r >= R
        if np.any(mask_ge_R):
            r_ge_R = r[mask_ge_R]
            v_sq[mask_ge_R] = V0 * (1.0 - (R / r_ge_R)**2)

        result = np.sqrt(v_sq)
        return result[0] if is_scalar else result

    # 3. Set up the radial Schrödinger equation for numerical solution (shooting method)

    def radial_ode_rhs(r, y, E):
        """
        Defines the system of first-order ODEs for the radial Schrödinger equation
        for an s-state (l=0). y is a vector [u, du/dr].
        """
        u, _ = y
        # The ODE is u''(r) = (V(r) - E) / (hbar^2/2m) * u(r)
        d2u_dr2 = (potential(r) - E) / hbar_sq_2m * u
        return [y[1], d2u_dr2]

    def shoot(E):
        """
        Solves the radial equation for a given energy E and returns the value
        of the radial wavefunction u(r) at a large distance r_max.
        The energy eigenvalues are the roots of this function (where u(r_max) -> 0).
        """
        r_max = 30.0    # Integration limit in nm, chosen to be far into the forbidden region
        r_start = 1e-9  # Start integration slightly away from r=0
        
        # For s-states (l=0), u(r) is proportional to r near the origin.
        # Initial conditions: u(r_start) = r_start, u'(r_start) = 1.0
        y0 = [r_start, 1.0]
        
        sol = solve_ivp(
            radial_ode_rhs,
            t_span=[r_start, r_max],
            y0=y0,
            args=(E,),
            method='RK45',
            dense_output=True
        )
        return sol.sol(r_max)[0]

    # 4. Find the first two energy levels (E1, E2) using a root-finder

    # Bound states must have energy E below the asymptotic potential V(inf) = sqrt(V0)
    # V(inf) is approx 3.87 eV.
    # From exploratory analysis, we can define brackets where the roots E1 and E2 lie.
    E1_bracket = [0.1, 1.5]  # Bracket for the ground state E1
    E2_bracket = [1.5, 3.5]  # Bracket for the first excited state E2

    E1 = brentq(shoot, E1_bracket[0], E1_bracket[1])
    E2 = brentq(shoot, E2_bracket[0], E2_bracket[1])
    
    # 5. Calculate the energy difference and print the final result

    delta_E = E2 - E1

    print("Calculation Results:")
    print("The first energy level (ground state) E1 = {:.4f} eV".format(E1))
    print("The second energy level (first excited state) E2 = {:.4f} eV".format(E2))
    print("\nThe energy difference between the levels is:")
    print("ΔE = E2 - E1 = {:.4f} eV - {:.4f} eV = {:.4f} eV".format(E2, E1, delta_E))
    
    return delta_E

if __name__ == '__main__':
    final_answer = solve_quantum_well_energy_difference()
    # The final value is printed within the function itself, as requested.
    # To conform to the output format, we also capture the return value.
    # print(f"<<<{final_answer:.4f}>>>")
