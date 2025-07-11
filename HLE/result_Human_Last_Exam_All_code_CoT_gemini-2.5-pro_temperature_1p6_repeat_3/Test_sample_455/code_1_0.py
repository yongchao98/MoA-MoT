import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize

# --- Constants and Parameters ---
# We use units of eV and nm.
# (hbar*c)^2 / (2*m_e*c^2)
HBARC = 197.32698  # eV nm
M_E_C2 = 511000   # eV (electron mass)
K_CONST = (HBARC**2) / (2 * M_E_C2) # This constant is approximately 0.038089 eV nm^2

# Parameters from the problem statement
V0_PARAM = 15.0   # Interpreted as a numerical value in the formula, yielding a potential in eV
R_PARAM = 3.0     # The radius in nm

# --- Numerical settings ---
R_MIN = 1e-5      # Start integration near r=0 to avoid singularities
R_MAX = 30.0      # A sufficiently large distance for the wavefunction to decay

def potential_eV(r):
    """
    Calculates the potential V(r) in eV, for r in nm.
    This function implements the interpretation of the problem statement's formula.
    """
    if r <= 0:
        return np.inf
    
    try:
        if r < R_PARAM:
            # Region r < R: V^2(r) = V0 + W(e^(r-R))
            v_sq = V0_PARAM + lambertw(np.exp(r - R_PARAM), k=0).real
        else: # r >= R_PARAM
            # Region r >= R: V^2(r) = V0 * (1 - (r/R)^-2)
            v_sq = V0_PARAM * (1.0 - (R_PARAM / r)**2)
        
        if v_sq < 0:
             return np.inf # Unphysical, treat as an infinite wall
        return np.sqrt(v_sq)

    except (ValueError, TypeError):
        return np.inf

def v_eff_eV(r, l):
    """Calculates the effective potential V_eff = V(r) + centrifugal term."""
    if r <= 0:
        return np.inf
    potential_val = potential_eV(r)
    centrifugal_term = K_CONST * l * (l + 1) / r**2
    return potential_val + centrifugal_term

def ode_system(r, y, E, l):
    """
    Defines the system of first-order ODEs for the shooting method.
    u'' = (V_eff - E) / K_CONST * u
    """
    u, du_dr = y
    potential_val = v_eff_eV(r, l)
    d2u_dr2 = (potential_val - E) / K_CONST * u
    return [du_dr, d2u_dr2]

def wave_function_at_rmax(E, l):
    """
    Solves the ODE for a given energy E and returns the wavefunction's value at R_MAX.
    This function's root is an energy eigenvalue.
    """
    u0 = R_MIN**(l + 1)
    du0_dr = (l + 1) * R_MIN**l
    y0 = [u0, du0_dr]
    
    sol = solve_ivp(
        fun=ode_system, t_span=[R_MIN, R_MAX], y0=y0,
        args=(E, l), dense_output=True, method='RK45'
    )
    return sol.sol(R_MAX)[0]

# The potential asymptotes to V(inf) = sqrt(V0)
V_inf = np.sqrt(V0_PARAM)

try:
    # Find E1 (ground state, n=1, l=0). Its energy must be between V_min=0 and V_inf.
    E1 = brentq(wave_function_at_rmax, 0.01, V_inf * 0.999, args=(0,))

    # Find the energy of the next l=0 state (n=2, l=0)
    E_2s = np.inf
    try:
        E_2s = brentq(wave_function_at_rmax, E1 + 1e-4, V_inf * 0.999, args=(0,))
    except ValueError:
        pass # This state may not exist

    # Find the energy of the lowest l=1 state (n=1, l=1)
    # The effective potential minimum for l=1 is > 0. Find it for a good search bracket.
    res = minimize(lambda r: v_eff_eV(r, 1), x0=R_PARAM, bounds=[(R_MIN, R_MAX)])
    min_v_eff_l1 = res.fun
    E_1p = np.inf
    try:
        # Search for the eigenvalue above the potential minimum.
        E_1p = brentq(wave_function_at_rmax, min_v_eff_l1 + 1e-4, V_inf * 0.999, args=(1,))
    except ValueError:
        pass # This state may not exist
        
    # The first excited state E2 is the minimum of the two possibilities
    E2 = min(E_2s, E_1p)

    if E2 == np.inf:
        print("Could not find a valid first excited state.")
    else:
        # Calculate the energy difference and print the final equation
        delta_E = E2 - E1
        print(f"{E2:.5f} - {E1:.5f} = {delta_E:.5f}")
        print(f"<<<{delta_E:.5f}>>>")

except ValueError:
    print("Could not find a bound ground state. The potential may not support any bound states with the given parameters.")
