import numpy as np
from scipy import constants
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_quantum_well_problem():
    """
    Solves for the energy levels of a particle in a 3D potential well
    and calculates the difference between the first two levels using the shooting method.
    """
    # 1. Constants and Parameters (in SI units)
    m = constants.m_e
    hbar = constants.hbar
    e = constants.e
    V0_eV = 15.0
    R_nm = 3.0
    V0 = V0_eV * e
    R = R_nm * 1e-9

    # 2. Potential Energy Function V(r)
    def potential(r):
        """Calculates V(r) from the given V^2(r) function."""
        # Ensure r is a numpy array for vectorized operations
        r = np.asarray(r)
        V_sq = np.zeros_like(r, dtype=float)

        # Region 1: 0 <= r < R
        mask_lt_R = (r >= 0) & (r < R)
        if np.any(mask_lt_R):
            r_lt_R = r[mask_lt_R]
            # Use .real to ensure the output is a real number
            V_sq[mask_lt_R] = V0 + lambertw(np.exp(r_lt_R - R)).real

        # Region 2: r >= R
        mask_ge_R = r >= R
        if np.any(mask_ge_R):
            r_ge_R = r[mask_ge_R]
            V_sq[mask_ge_R] = V0 * (1 - (R / r_ge_R)**2)
        
        # V(r) is the square root of V^2(r)
        return np.sqrt(V_sq)

    # 3. Radial ODE system
    def radial_ode(r, y, E, l):
        """Defines the system of ODEs for the radial Schrödinger equation."""
        u, u_prime = y
        # The centrifugal term l(l+1)/r^2
        # An exception is added for r=0, though the integration starts at r_min > 0.
        if r == 0:
            centrifugal_term = 0
        else:
            centrifugal_term = l * (l + 1) / r**2
            
        d_u_prime_dr = ((2 * m / hbar**2) * (potential(r) - E) + centrifugal_term) * u
        return [u_prime, d_u_prime_dr]

    # 4. Shooting method function
    def wave_func_at_rmax(E, l, r_span, r_min):
        """
        Solves the ODE for a given energy E and returns the value of u(r_max).
        The roots of this function in E are the energy eigenvalues.
        """
        # Initial conditions at r_min > 0, consistent with u(r) ~ r^(l+1) near the origin.
        y0 = [r_min**(l + 1), (l + 1) * r_min**l]
        
        # Solve the initial value problem with specified tolerances
        sol = solve_ivp(radial_ode, r_span, y0, args=(E, l), method='RK45', atol=1e-8, rtol=1e-8)
        
        # Return the value of the wavefunction at the last point (r_max)
        return sol.y[0, -1]

    # 5. Eigenvalue finding function
    def find_eigenvalues(l, num_levels_to_find, E_range_J, r_span, r_min):
        """
        Scans an energy range to find brackets for eigenvalues, then uses a
        root-finder to get the precise values.
        """
        E_min_J, E_max_J = E_range_J
        eigenvalues = []
        
        # Create a grid of energies to scan for sign changes
        test_energies = np.linspace(E_min_J, E_max_J, 700)
        u_at_rmax_values = np.array([wave_func_at_rmax(E, l, r_span, r_min) for E in test_energies])

        # Find indices where the function u(r_max) crosses the axis (a sign change)
        sign_change_indices = np.where(np.diff(np.sign(u_at_rmax_values)))[0]

        for i in sign_change_indices:
            if len(eigenvalues) >= num_levels_to_find:
                break
            
            E_low_bracket = test_energies[i]
            E_high_bracket = test_energies[i+1]
            
            try:
                # Find the root (eigenvalue) within the identified bracket
                eigenvalue = brentq(wave_func_at_rmax, E_low_bracket, E_high_bracket, args=(l, r_span, r_min))
                eigenvalues.append(eigenvalue)
            except ValueError:
                continue # brentq fails if signs are not opposite, skip bracket.
                
        return np.array(eigenvalues)

    # --- Main execution of the plan ---
    # Integration range settings
    r_min = 1e-15  # Start integration slightly away from the origin
    r_max = 20 * R # A sufficiently large radius where wavefunction should have decayed
    r_span = [r_min, r_max]

    # Energy search range for bound states (E > V_min, E < V_asymptote)
    # V_min = V(R) = 0. V_asymptote = V(inf) = sqrt(V0).
    # Search for energies slightly above the potential minimum and below the asymptote.
    E_min_J = 0.01 * e
    E_max_J = np.sqrt(V0) - (0.01 * e)

    # 6. Find necessary energy levels
    # Find the two lowest energy levels for l=0 (s-states)
    energies_l0 = find_eigenvalues(l=0, num_levels_to_find=2, E_range_J=(E_min_J, E_max_J), r_span=r_span, r_min=r_min)
    if len(energies_l0) < 2:
        print("Error: Could not find two energy levels for l=0. The potential well may be too shallow or the search range incorrect.")
        return

    E_1_0_J = energies_l0[0]
    E_2_0_J = energies_l0[1]

    # Find the lowest energy level for l=1 (p-state)
    energies_l1 = find_eigenvalues(l=1, num_levels_to_find=1, E_range_J=(E_min_J, E_max_J), r_span=r_span, r_min=r_min)
    if len(energies_l1) < 1:
        E_1_1_J = np.inf # If no l=1 state found, set its energy to infinity
    else:
        E_1_1_J = energies_l1[0]
    
    # 7. Determine the overall ground state energy, E1
    E1_J = E_1_0_J

    # 8. Determine the overall first excited state energy, E2
    E2_J = min(E_2_0_J, E_1_1_J)
    
    # 9. Calculate the energy difference
    delta_E_J = E2_J - E1_J

    # 10. Convert to eV and print the result
    E1_eV = E1_J / e
    E2_eV = E2_J / e
    delta_E_eV = delta_E_J / e

    print("The energy difference between the first and second energy levels is calculated as follows:")
    print(f"ΔE = E₂ - E₁ = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")

# Execute the solver function
solve_quantum_well_problem()
<<<0.4287>>>