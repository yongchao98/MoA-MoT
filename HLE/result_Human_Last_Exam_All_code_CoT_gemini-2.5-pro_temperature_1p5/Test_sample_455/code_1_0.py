import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# 1. Define physical constants and parameters in consistent units (eV, nm).
V0 = 15.0  # eV
R = 3.0    # nm
# We need the constant hbar^2 / (2m)
HBAR_C = 197.32698  # eV*nm
M_E_C2 = 510998.95 # eV (mass of electron in energy units)
HBAR2_OVER_2M = (HBAR_C**2) / (2 * M_E_C2) # This is ~0.0380998 eV*nm^2

# 2. Implement the potential energy function V(r) from V^2(r).
def potential(r_nm):
    """
    Calculates the potential energy V(r) in eV for a given radius r in nm.
    """
    # The potential is defined piecewise.
    if r_nm < R:
        # V^2(r) = V0 + W(exp(r-R)). For r < R, the argument to exp is negative.
        # W(exp(x)) for x<0 is real and can be computed by lambertw.
        potential_squared = V0 + lambertw(np.exp(r_nm - R)).real
    else:  # r_nm >= R
        # V^2(r) = V0 * (1 - (r/R)^-2). This simplifies to V0*(1 - R^2/r^2).
        # This term is always non-negative for r>=R.
        potential_squared = V0 * (1.0 - (r_nm / R)**(-2.0))
            
    return np.sqrt(potential_squared)

# 3. Define the ODE system for the radial Schrödinger equation.
def schrodinger_ode(r, y, E, l):
    """
    Defines the system of first-order ODEs for u(r).
    y is the vector [u(r), u'(r)]. E is energy, l is angular momentum quantum number.
    """
    u, du_dr = y
    
    # The centrifugal term is l(l+1) hbar^2 / (2mr^2)
    # Handle r=0 to avoid division by zero, though we start integration at r>0.
    if r == 0:
        centrifugal_term = 0.0
    else:
        centrifugal_term = l * (l + 1) * HBAR2_OVER_2M / (r**2)
    
    effective_potential = potential(r) + centrifugal_term
    d2u_dr2 = - (E - effective_potential) * u / HBAR2_OVER_2M
    return [du_dr, d2u_dr2]

# 4. Implement the shooting method logic.
def shoot_wavefunction(E, l):
    """
    Solves the ODE for a given energy E and quantum number l.
    Returns the value of u(r) at a large radius r_end.
    Eigenvalues are the roots of this function.
    """
    r_start = 1e-6  # Start integration slightly away from the origin
    r_end = 8 * R   # Integrate out to a large radius in the forbidden region

    # Initial conditions for u(r) ~ r^(l+1) near r=0
    y_start = [r_start**(l+1), (l+1)*r_start**l]

    sol = solve_ivp(
        fun=lambda r, y: schrodinger_ode(r, y, E, l),
        t_span=[r_start, r_end],
        y0=y_start,
        t_eval=[r_end] # Only need the final value
    )
    
    return sol.y[0, -1]

# 5. Find the first few energy eigenvalues for a given l.
def find_energy_levels(l, num_levels_to_find):
    """
    Scans for the first few energy eigenvalues for a given l by finding roots of shoot_wavefunction.
    """
    levels = []
    # Bound states must have energy below the potential at infinity.
    V_inf = np.sqrt(V0)
    
    # Start searching for energy just above the minimum of the effective potential.
    if l == 0:
        E_min_search = 0.01
    else:
        # For l>0, the centrifugal barrier creates a minimum > 0.
        # This is a rough estimate; we start searching from a small positive energy.
        E_min_search = 0.1 

    e_scan = np.linspace(E_min_search, V_inf * 0.99, 1000)
    u_end_values = [shoot_wavefunction(e, l) for e in e_scan]

    # Find roots by locating sign changes
    for i in range(len(e_scan) - 1):
        if len(levels) >= num_levels_to_find:
            break
        if np.sign(u_end_values[i]) != np.sign(u_end_values[i+1]):
            e_low, e_high = e_scan[i], e_scan[i+1]
            try:
                energy = brentq(lambda e: shoot_wavefunction(e, l), e_low, e_high)
                levels.append(energy)
            except (ValueError, RuntimeError):
                # brentq can fail if the function is strange, just continue.
                continue
    return levels

# 6. Main execution: Calculate the required energy levels and their difference.
if __name__ == "__main__":
    # Find the lowest two s-wave (l=0) states.
    s_levels = find_energy_levels(l=0, num_levels_to_find=2)
    E_1s = s_levels[0] if len(s_levels) > 0 else float('inf')
    E_2s = s_levels[1] if len(s_levels) > 1 else float('inf')
    
    # Find the lowest p-wave (l=1) state.
    p_levels = find_energy_levels(l=1, num_levels_to_find=1)
    E_1p = p_levels[0] if len(p_levels) > 0 else float('inf')
    
    # The first energy level (ground state) is E_1s.
    E1 = E_1s
    
    # The second energy level is the minimum of the next available states.
    E2 = min(E_2s, E_1p)
    
    if E1 == float('inf') or E2 == float('inf'):
        print("Could not find at least two bound states for this potential.")
    else:
        # Calculate the energy difference.
        delta_E = E2 - E1
        
        # Print the results as requested.
        print("The first energy level (ground state) E1 is {:.4f} eV.".format(E1))
        print("The second energy level (first excited state) E2 is {:.4f} eV.".format(E2))
        print("The energy difference is Delta E = E2 - E1")
        print("Delta E = {:.4f} eV - {:.4f} eV = {:.4f} eV".format(E2, E1, delta_E))
        # Final answer token for parsing.
        # print(f"<<<{delta_E:.4f}>>>")