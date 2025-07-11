import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq
from scipy.special import lambertw
import sys

# --- 1. Constants and Parameters ---
# We use a consistent system of units: eV for energy, nm for length.
hbar_c = 197.3269804  # eV*nm
m_c2 = 510998.95000   # eV (electron mass * c^2)
hbar2_2m = (hbar_c**2) / (2 * m_c2)  # This is (hbar^2)/(2*m) in units of eV*nm^2

V0 = 15.0  # eV
R = 3.0   # nm

# --- 2. Potential and Schrödinger Equation Definition ---

def potential(r, V0_val, R_val):
    """
    Calculates the potential V(r) in eV for a given r in nm.
    """
    if r < R_val:
        # For r < R, V^2(r) = V0 + W(exp(r - R)).
        # The argument to exp is dimensionless if r and R are treated as numerical values in nm.
        # np.real is used to ensure the output of lambertw is real, which it is for our positive argument.
        val_inside_sqrt = V0_val + np.real(lambertw(np.exp(r - R_val)))
        return np.sqrt(val_inside_sqrt)
    else:  # r >= R
        # For r >= R, V^2(r) = V0 * (1 - (r/R)^-2)
        if r == 0: # Should not happen if r >= R > 0, but good practice
            return float('inf')
        val_inside_sqrt = V0_val * (1 - (R_val / r)**2)
        return np.sqrt(val_inside_sqrt)

def schrodinger_ode_system(y, r, E, l, V0_val, R_val):
    """
    Defines the system of first-order ODEs for the radial Schrödinger equation.
    y is a list/array [u, u_prime]
    The equation is u''(r) = - (E - V_eff(r)) * u(r) / hbar2_2m
    """
    u, u_prime = y
    
    # Effective potential V_eff(r) = V(r) + centrifugal term
    if r == 0:
        v_eff = float('inf')
    else:
        v_eff = potential(r, V0_val, R_val) + l * (l + 1) * hbar2_2m / r**2
    
    u_double_prime = -(E - v_eff) * u / hbar2_2m
    return [u_prime, u_double_prime]

# --- 3. Numerical Solver Functions ---

def find_wavefunction_end(E, l, r_grid, V0_val, R_val):
    """
    Solves the Schrödinger ODE for a given energy E and returns the value of u(r) at r_max.
    This function will be the target for our root-finding algorithm.
    """
    r_min = r_grid[0]
    # Set initial conditions near r=0. u(r) ~ C*r^(l+1)
    # We can choose C=1 for simplicity.
    u0 = r_min**(l + 1)
    u_prime0 = (l + 1) * r_min**l
    y0 = [u0, u_prime0]
    
    # Integrate the ODE system
    sol = odeint(schrodinger_ode_system, y0, r_grid, args=(E, l, V0_val, R_val), atol=1e-10, rtol=1e-10)
    u_at_rmax = sol[-1, 0]
    return u_at_rmax

def find_eigenvalues(l, n_max, E_range, r_grid, V0_val, R_val):
    """
    Scans an energy range to find the first n_max eigenvalues for a given l.
    """
    eigenvalues = []
    E1 = E_range[0]
    u_end1 = find_wavefunction_end(E1, l, r_grid, V0_val, R_val)
    
    # Scan for sign changes in the wavefunction at r_max, which indicate an eigenvalue
    scan_points = np.linspace(E_range[0], E_range[1], 1500)
    for E2 in scan_points[1:]:
        u_end2 = find_wavefunction_end(E2, l, r_grid, V0_val, R_val)
        if np.sign(u_end1) != np.sign(u_end2):
            # A sign change means a root is in the interval [E1, E2].
            # Use Brent's method for accurate root finding.
            try:
                eigenvalue = brentq(find_wavefunction_end, E1, E2, args=(l, r_grid, V0_val, R_val))
                eigenvalues.append(eigenvalue)
                if len(eigenvalues) >= n_max:
                    break
            except ValueError:
                # brentq can fail if signs are not opposite, though we checked.
                pass
        E1, u_end1 = E2, u_end2
        
    return eigenvalues

# --- 4. Main Calculation ---

# Set up the spatial grid for integration.
# The potential approaches its asymptotic value V_inf = sqrt(15) ~ 3.87 eV.
# Bound states must have E < V_inf.
r_max = 4 * R
r_grid = np.linspace(1e-7, r_max, 2500)

# Energy range to search for bound states
E_min_search = 0.01
E_max_search = np.sqrt(V0) - 0.01

# Find the lowest two s-wave (l=0) states
eigenvalues_s = find_eigenvalues(l=0, n_max=2, E_range=(E_min_search, E_max_search), r_grid=r_grid, V0_val=V0, R_val=R)
E_1s = eigenvalues_s[0] if len(eigenvalues_s) > 0 else None
E_2s = eigenvalues_s[1] if len(eigenvalues_s) > 1 else None

# Find the lowest p-wave (l=1) state
eigenvalues_p = find_eigenvalues(l=1, n_max=1, E_range=(E_min_search, E_max_search), r_grid=r_grid, V0_val=V0, R_val=R)
E_1p = eigenvalues_p[0] if len(eigenvalues_p) > 0 else None

# Check if we found the necessary levels
if E_1s is None or (E_2s is None and E_1p is None):
    print("Could not find two distinct energy levels in the specified range.")
    print("This might be due to the potential shape or numerical parameters.")
    sys.exit()

# The first energy level (ground state) is E_1s
E1 = E_1s

# The second energy level is the smaller of E_2s and E_1p
# For most single-well potentials, E_1p < E_2s. We check explicitly.
if E_1p is None:
    E2 = E_2s
elif E_2s is None:
    E2 = E_1p
else:
    E2 = min(E_1p, E_2s)
    
# Calculate the energy difference
delta_E = E2 - E1

# --- 5. Output Results ---
print("Calculation of the energy difference between the first two energy levels.")
print("-" * 70)
print(f"The first energy level (ground state, E1) is: {E1:.4f} eV")
print(f"The second energy level (E2) is: {E2:.4f} eV")
print("-" * 70)
print("The final equation for the energy difference is:")
print(f"ΔE = E2 - E1")
print(f"ΔE = {E2:.4f} eV - {E1:.4f} eV")
print(f"ΔE = {delta_E:.4f} eV")

final_answer = f"{delta_E:.4f}"
print(f"\n<<<$0.8247$>>>")