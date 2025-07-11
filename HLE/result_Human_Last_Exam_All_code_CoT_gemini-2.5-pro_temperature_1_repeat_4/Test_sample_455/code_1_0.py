import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# This script calculates the first two energy levels of a quantum particle
# in a given 3D potential well and finds their difference.

# 1. Define physical constants and problem parameters.
# All energy units are in eV, and distance units are in nm.
V0 = 15.0  # eV
R = 3.0    # nm
hbar_c = 197.327  # eV·nm
m_c2 = 511000     # electron rest mass energy in eV
# The constant term in the Schrödinger equation: ħ²/2m = (ħc)²/2mc²
C = hbar_c**2 / (2 * m_c2)  # eV·nm²

# 2. Define the potential energy function U(r).
# The problem uses the label V²(r) for the potential energy function.
def potential_energy(r):
    """
    Calculates the potential energy U(r) in eV.
    """
    if r < R:
        # For r < R, U(r) = V₀ + W(exp(r - R))
        # The principal branch of Lambert W function is real for positive arguments.
        return V0 + np.real(lambertw(np.exp(r - R)))
    else:
        # For r ≥ R, U(r) = V₀ * (1 - (R/r)²)
        return V0 * (1 - (R / r)**2)

# 3. Set up the ODE system for the shooting method.
# The radial Schrödinger equation is u''(r) = - (E - U(r)) / C * u(r)
def ode_system(r, y, E):
    """
    Defines the system of first-order ODEs.
    y = [u, u'] -> dy/dr = [u', u'']
    """
    u, u_prime = y
    u_double_prime = -(E - potential_energy(r)) / C * u
    return [u_prime, u_double_prime]

# 4. Define the function whose roots are the energy eigenvalues.
def find_wave_function_end_value(E):
    """
    Solves the ODE for a given energy E and returns the value of the
    wavefunction u at a large radius, r_max.
    """
    r_min = 1e-6    # Start integration slightly away from r=0
    r_max = 6 * R   # Integrate to a sufficiently large radius (18 nm)
    
    # Boundary conditions at r=0: u(0)=0. We set u'(0)=1.
    y0 = [0.0, 1.0]
    
    # Numerically solve the initial value problem
    sol = solve_ivp(
        fun=lambda r, y: ode_system(r, y, E),
        t_span=[r_min, r_max],
        y0=y0,
        dense_output=True,
        method='RK45'
    )
    
    # Return the value of u at r_max. Eigenvalues are roots of this function.
    return sol.y[0, -1]

# 5. Find the energy eigenvalues (E1, E2) by finding roots.
# Bound states must have energy 0 < E < V₀. We scan this range.
energy_range = np.linspace(0.1, V0 - 0.1, 500)
wave_func_values = [find_wave_function_end_value(e) for e in energy_range]

eigenvalues = []
# Find intervals where the function crosses zero, then find the root in that interval.
for i in range(len(energy_range) - 1):
    if np.sign(wave_func_values[i]) != np.sign(wave_func_values[i+1]):
        e_low = energy_range[i]
        e_high = energy_range[i+1]
        try:
            root = brentq(find_wave_function_end_value, e_low, e_high)
            eigenvalues.append(root)
        except (ValueError, RuntimeError):
            pass # Ignore cases where root finding fails

if len(eigenvalues) < 2:
    print("Could not find two energy levels in the specified range.")
else:
    E1 = eigenvalues[0]
    E2 = eigenvalues[1]
    delta_E = E2 - E1

    # 6. Print the final results in the required format.
    print(f"The first energy level (ground state), E1 = {E1:.4f} eV")
    print(f"The second energy level (first excited state), E2 = {E2:.4f} eV")
    print(f"The energy difference is ΔE = E2 - E1")
    print(f"ΔE = {E2:.4f} eV - {E1:.4f} eV = {delta_E:.4f} eV")
    
    print(f"\n<<<The final calculated energy difference is {delta_E:.4f} eV.>>>")
