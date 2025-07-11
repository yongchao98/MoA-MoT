import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq

def solve_quantum_well():
    """
    Calculates the energy difference between the first two energy levels
    of a particle in the given potential well.
    """
    # --- Constants in eV and nm ---
    # hbar^2 / (2 * m_e)
    HBAR2_2M = 0.0380998  # eV * nm^2
    V0 = 15.0  # eV
    R = 3.0  # nm

    # --- Potential Energy Function ---
    def V_eff(r, l):
        """
        Calculates the effective potential for r >= R.
        The region r < R is treated as an infinite wall.
        """
        if r <= R:
            return np.inf
        # Centrifugal term
        centrifugal_term = HBAR2_2M * l * (l + 1) / (r**2)
        # Potential term from the problem
        potential_term = np.sqrt(V0 * (1 - (R / r)**2))
        return centrifugal_term + potential_term

    # --- Schrödinger ODE System ---
    def schrodinger_ode(y, r, E, l):
        """
        Defines the system of first-order ODEs for the radial Schrödinger equation.
        y = [u, u'], where u is the radial wavefunction.
        Returns [u', u''].
        """
        u, du_dr = y
        u_double_prime = (V_eff(r, l) - E) / HBAR2_2M * u
        return [du_dr, u_double_prime]

    # --- Shooting Method Objective Function ---
    def objective_function(E, l, r_max, dr):
        """
        Solves the ODE for a given energy E and returns the value of the
        wavefunction u at a large radius r_max.
        The energy eigenvalues are the roots of this function.
        """
        # Start integration slightly away from R to avoid singularity in V'(R)
        r_start = R + dr
        # Initial conditions: u(R)=0, u'(R)=1 (arbitrary slope)
        y0 = [dr, 1.0]
        # Create the array of r values for the integration
        r_span = np.arange(r_start, r_max, dr)
        # Solve the ODE
        solution = odeint(schrodinger_ode, y0, r_span, args=(E, l), tfirst=True)
        # Return the wavefunction's value at the endpoint
        return solution[-1, 0]

    # --- Numerical Parameters ---
    r_max = 10.0  # Maximum radius for integration (nm)
    dr = 0.002    # Step size for integration (nm)
    E_min = 1e-6  # Minimum energy to search from (eV)
    E_max = np.sqrt(V0) - 1e-6 # Max energy is the asymptotic potential (eV)

    # --- Find Energy Eigenvalues ---

    # Find E_1,0 (ground state, n=1, l=0)
    E_1_0 = brentq(objective_function, E_min, E_max, args=(0, r_max, dr))

    # Find E_1,1 (first l=1 state, n=1, l=1)
    # This state must have higher energy than E_1_0
    E_1_1 = brentq(objective_function, E_1_0, E_max, args=(1, r_max, dr))
    
    # Find E_2,0 (first excited state for l=0, n=2, l=0)
    # We need to find the second root for l=0. We search for the next sign change.
    E_search_start = E_1_0 + 0.01
    e_low = E_search_start
    e_step = 0.01
    sign_low = np.sign(objective_function(e_low, 0, r_max, dr))
    E_2_0 = None
    while e_low < E_max:
        e_high = e_low + e_step
        if e_high > E_max:
            break
        sign_high = np.sign(objective_function(e_high, 0, r_max, dr))
        if sign_high != sign_low:
            E_2_0 = brentq(objective_function, e_low, e_high, args=(0, r_max, dr))
            break
        e_low = e_high
        sign_low = sign_high

    # --- Determine E1 and E2 and Calculate Delta E ---
    E1 = E_1_0
    
    # The second energy level is the minimum of the next available states
    if E_2_0 is not None:
        E2 = min(E_1_1, E_2_0)
    else: # Should not happen in this potential
        E2 = E_1_1

    delta_E = E2 - E1

    # --- Print the Final Results ---
    print(f"Ground state energy E1 (n=1, l=0): {E1:.4f} eV")
    print(f"First excited state candidates:")
    print(f"  - E(n=1, l=1): {E_1_1:.4f} eV")
    if E_2_0 is not None:
        print(f"  - E(n=2, l=0): {E_2_0:.4f} eV")
    print(f"Second energy level E2 is: {E2:.4f} eV")
    print("\nCalculation of the energy difference:")
    print(f"Delta E = E2 - E1 = {E2:.4f} - {E1:.4f} = {delta_E:.4f} eV")
    print(f"\n<<<{delta_E:.4f}>>>")

solve_quantum_well()
