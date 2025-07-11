import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_quantum_well():
    """
    Solves for the energy levels of a particle in the given 3D potential well.
    """
    # --- 1. Define constants and parameters in consistent units (eV, nm) ---
    # Mass of the particle (electron) in eV/c^2
    M_EV_C2 = 0.511e6  
    # Planck constant-c in eV*nm
    HBARC = 197.327
    # ħ^2 / 2m in eV * nm^2
    HBAR_SQ_2M = (HBARC**2) / (2 * M_EV_C2)

    # Potential parameters
    V0 = 15.0  # eV
    R = 3.0    # nm

    # --- 2. Define the potential energy function V(r) ---
    def potential(r):
        """
        Calculates the potential V(r) in eV for r in nm.
        Handles the piecewise definition of the potential.
        """
        if r <= 1e-9:  # Avoid r=0 singularity and handle very small r
            return np.inf
        
        try:
            if r < R:
                # For r < R, V^2 = V0 + W(exp(r - R))
                # The principal branch of Lambert W is used (k=0).
                # np.real is a safeguard; for positive arguments, the result is real.
                arg_w = np.exp(r - R)
                v_sq = V0 + np.real(lambertw(arg_w, k=0))
            else:  # r >= R
                # For r >= R, V^2 = V0 * (1 - (r/R)^-2)
                v_sq = V0 * (1.0 - (R / r)**2)
            
            # The potential is the square root of V^2
            # A negative v_sq would lead to imaginary potential (not the case here)
            return np.sqrt(v_sq)
        except (ValueError, TypeError):
            return np.inf

    # --- 3. Define the ODE for the numerical solver ---
    def radial_ode(r, y, E, l):
        """
        Represents the radial Schrödinger equation as a system of two 1st order ODEs.
        y is a vector [u, du/dr].
        Returns [du/dr, d^2u/dr^2].
        """
        u, du_dr = y
        
        # Effective potential V_eff(r) = V(r) + centrifugal term
        V_eff = potential(r) + HBAR_SQ_2M * l * (l + 1) / (r**2)
        
        # d^2u/dr^2 = (2m/ħ^2) * (V_eff - E) * u
        d2u_dr2 = (V_eff - E) / HBAR_SQ_2M * u
        
        return [du_dr, d2u_dr2]

    # --- 4. Function to find energy levels using the shooting method ---
    def find_energy_level(principal_n, l):
        """
        Finds the n-th energy eigenvalue for a given angular momentum l.
        principal_n=1 is the lowest energy state for that l.
        """
        # Integration range [r_min, r_max]
        r_min = 1e-5  # Start integration slightly away from r=0
        r_max = 8 * R   # Integrate to a point where the wavefunction should be zero

        # Boundary conditions at r_min, based on u(r) ~ r^(l+1) near the origin
        u_min = r_min**(l + 1)
        du_min = (l + 1) * r_min**l
        initial_conditions = [u_min, du_min]
        
        # The objective function whose roots are the energy eigenvalues.
        # It's the value of the wavefunction at r_max.
        def objective_function(E):
            sol = solve_ivp(
                radial_ode,
                [r_min, r_max],
                initial_conditions,
                args=(E, l),
                method='RK45'
            )
            return sol.y[0, -1]

        # Search for energy intervals where the wavefunction sign flips at r_max
        E_min_search = 1e-6
        E_max_search = np.sqrt(V0) * 0.999 # Max energy for a bound state
        
        energy_scan = np.linspace(E_min_search, E_max_search, 400)
        u_values = [objective_function(E) for E in energy_scan]
        
        roots = []
        for i in range(len(u_values) - 1):
            if np.sign(u_values[i]) != np.sign(u_values[i+1]):
                E_low, E_high = energy_scan[i], energy_scan[i+1]
                try:
                    root = brentq(objective_function, E_low, E_high)
                    roots.append(root)
                except ValueError:
                    pass
        
        if principal_n > len(roots):
            return None # The requested energy level was not found
        
        return roots[principal_n - 1]

    # --- 5. Calculate the required energy levels and the difference ---
    print("Calculating energy levels...")
    
    # E1 is the ground state (n=1, l=0)
    E1 = find_energy_level(principal_n=1, l=0)
    
    # For E2, we find the next lowest states for l=0 and l=1
    E2_s_state = find_energy_level(principal_n=2, l=0) # second s-state
    E2_p_state = find_energy_level(principal_n=1, l=1) # first p-state
    
    print("\n--- Results ---")
    if E1 is None:
        print("Could not find the ground state E1.")
        return

    print(f"Ground state energy (E1), for n=1, l=0: E_10 = {E1:.4f} eV")

    valid_excited_states = []
    if E2_s_state is not None:
        valid_excited_states.append(E2_s_state)
        print(f"First excited s-state energy, for n=2, l=0: E_20 = {E2_s_state:.4f} eV")
    if E2_p_state is not None:
        valid_excited_states.append(E2_p_state)
        print(f"Lowest p-state energy, for n=1, l=1: E_11 = {E2_p_state:.4f} eV")

    if not valid_excited_states:
        print("Could not find any excited states.")
        return
        
    # E2 is the minimum of the found excited states
    E2 = min(valid_excited_states)
    
    Delta_E = E2 - E1

    print("\n--- Final Calculation ---")
    print(f"The first energy level E1 is {E1:.4f} eV.")
    print(f"The second energy level E2 is the minimum of the calculated excited states, so E2 = {E2:.4f} eV.")
    print(f"The energy difference is ΔE = E2 - E1 = {E2:.4f} - {E1:.4f} = {Delta_E:.4f} eV.")
    
# Execute the solver function
solve_quantum_well()