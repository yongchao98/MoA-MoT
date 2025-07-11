import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_quantum_well_problem():
    """
    Calculates the energy difference between the first and second energy levels
    of a particle in a 3D potential well by numerically solving the Schrödinger equation.
    """

    # --- Step 1: Define Constants and the Potential ---

    # Physical constants in units of eV and nm
    V0 = 15.0  # eV
    R = 3.0    # nm
    HBARC = 197.3  # eV * nm
    MEC2 = 0.511e6 # electron rest mass energy in eV
    HBAR2_2M = (HBARC**2) / (2 * MEC2)  # Approx. 0.0381 eV * nm^2

    # The potential function V(r) is the square root of the given V^2(r)
    def potential(r):
        """Calculates potential energy V(r) in eV."""
        if r < R:
            # For r<R, argument to lambertw is exp(r-R) which is > 0.
            # The principal branch W_0 is real, so np.real is a safeguard.
            potential_squared = V0 + np.real(lambertw(np.exp(r - R)))
        else:  # r >= R
            potential_squared = V0 * (1 - (R / r)**2)
        
        # Ensure potential_squared is not negative due to float precision
        return np.sqrt(max(0, potential_squared))

    # --- Step 2: Implement the Numerical Solver (Shooting Method) ---

    def get_log_derivative_mismatch(E, l):
        """
        Calculates the mismatch in the logarithmic derivative of the wavefunction
        at the matching radius r=R. Roots of this function are eigenvalues.
        """
        # Bound states must have energy less than the potential at infinity
        V_inf = np.sqrt(V0)
        if E >= V_inf:
            return np.inf  # Not a bound state

        # Define the ODE system y'' = f(r)*y, where y=[u, u']
        def ode_system(r, y):
            u, du_dr = y
            # Effective potential includes the centrifugal barrier
            V_eff = potential(r) + l * (l + 1) * HBAR2_2M / (r**2)
            d2u_dr2 = (V_eff - E) / HBAR2_2M * u
            return [du_dr, d2u_dr2]

        # Left integration: from near r=0 to r=R
        r_min = 1e-6  # A small radius to avoid r=0 singularity
        # Initial conditions for u(r) ~ r^(l+1) near the origin
        u_init_left = [r_min**(l + 1), (l + 1) * r_min**l]
        sol_left = solve_ivp(ode_system, [r_min, R], u_init_left, dense_output=True)
        u_left_R, du_left_R = sol_left.sol(R)
        if u_left_R == 0: return np.inf # Avoid division by zero
        log_deriv_left = du_left_R / u_left_R

        # Right integration: from r_max (large r) inwards to r=R
        r_max = R * 10
        # For large r, u(r) ~ exp(-kappa*r), where kappa determines the decay rate
        kappa = np.sqrt(max(0, (V_inf - E)) / HBAR2_2M)
        # Initial conditions for decaying solution at r_max (normalization is irrelevant)
        u_init_right = [np.exp(-kappa * r_max), -kappa * np.exp(-kappa * r_max)]
        sol_right = solve_ivp(ode_system, [r_max, R], u_init_right, dense_output=True)
        u_right_R, du_right_R = sol_right.sol(R)
        if u_right_R == 0: return np.inf # Avoid division by zero
        log_deriv_right = du_right_R / u_right_R
        
        # Return the difference
        return log_deriv_left - log_deriv_right

    def find_energy_levels(l, num_levels):
        """Finds the first 'num_levels' eigenvalues for a given l."""
        eigenvalues = []
        # Define the search range for energy. Min energy must be > V_eff minimum.
        # Max energy must be < V(inf).
        E_min = 0.01
        E_max = np.sqrt(V0) - 0.01
        
        # Scan the energy range to find brackets for the root solver
        energies = np.linspace(E_min, E_max, 500)
        mismatch_values = [get_log_derivative_mismatch(E, l) for E in energies]
        
        # Find roots where the sign of the mismatch function changes
        for i in range(len(mismatch_values) - 1):
            if np.sign(mismatch_values[i]) != np.sign(mismatch_values[i+1]):
                e1, e2 = energies[i], energies[i+1]
                try:
                    root = brentq(lambda E: get_log_derivative_mismatch(E, l), e1, e2)
                    eigenvalues.append(root)
                    if len(eigenvalues) == num_levels:
                        break
                except (ValueError, RuntimeError):
                    continue
        return eigenvalues

    # --- Step 3: Find Relevant Energy Levels ---
    
    # Find the two lowest s-wave (l=0) levels
    s_levels = find_energy_levels(l=0, num_levels=2)
    E_1s = s_levels[0]
    E_2s = s_levels[1]

    # Find the lowest p-wave (l=1) level
    p_levels = find_energy_levels(l=1, num_levels=1)
    E_1p = p_levels[0]
    
    # --- Step 4: Determine E1, E2, and Delta E ---

    # The ground state energy (E1) is the lowest s-wave state
    E1 = E_1s
    
    # The first excited state energy (E2) is the lower of the second s-wave and lowest p-wave states
    E2 = min(E_2s, E_1p)
    
    # Calculate the energy difference
    delta_E = E2 - E1

    # --- Step 5: Print the Final Result ---

    print(f"Ground state energy (1s): E1 = {E1:.6f} eV")
    if E2 == E_2s:
        print(f"First excited state energy (2s): E2 = {E2:.6f} eV")
    else:
        print(f"First excited state energy (1p): E2 = {E2:.6f} eV")
    
    print("\nFinal Calculation:")
    print(f"Energy difference: ΔE = E2 - E1 = {E2:.6f} eV - {E1:.6f} eV = {delta_E:.6f} eV")

    # Final answer in the required format
    print(f"\n<<<{delta_E:.6f}>>>")

# Execute the full calculation and print the result
solve_quantum_well_problem()