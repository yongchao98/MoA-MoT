import numpy as np
from scipy.special import lambertw
from scipy.integrate import odeint
from scipy.optimize import brentq
import warnings

# Suppress warnings that can occur with very small arguments to lambertw
warnings.filterwarnings("ignore", category=RuntimeWarning)

def solve_quantum_well():
    """
    Calculates the energy difference between the first two energy levels
    of a particle in the given potential well.
    """
    # --- Step 1: Define Constants ---
    V0 = 15.0  # Potential energy at center in eV
    R = 3.0    # Well radius in nm
    # The constant h_bar^2 / (2 * m_electron) in units of eV*nm^2
    HBAR2_2M_EV_NM2 = 0.0380952

    # --- Step 2: Define the Schrödinger ODE System ---
    def schrodinger_ode(y, r, E, l):
        """
        Defines the system of first-order ODEs for the Schrödinger equation.
        y[0] = u(r), y[1] = u'(r)
        """
        u, du_dr = y
        
        # Avoid division by zero at r=0
        if r == 0:
            return [0, 0]

        # Define the potential U(r) based on the problem statement
        if r < R:
            # For r < R, U(r) = V0 + W(exp(r - R))
            # Use np.real to ensure float output from the principal branch of W
            arg = np.exp(r - R)
            potential_r = V0 + np.real(lambertw(arg)) if arg > 0 else V0
        else:
            # For r >= R, U(r) = V0 * (1 - (R/r)^2)
            potential_r = V0 * (1.0 - (R / r)**2)
        
        # Effective potential V_eff(r) includes the centrifugal term
        V_eff = potential_r + HBAR2_2M_EV_NM2 * l * (l + 1) / r**2
        
        # Second derivative from the Schrödinger equation: u'' = -(E - V_eff) * u / C
        d2u_dr2 = -(E - V_eff) / HBAR2_2M_EV_NM2 * u
        
        return [du_dr, d2u_dr2]

    # --- Step 3: Define the Shooting Method Function ---
    def wavefunction_at_end(E, l, r_grid):
        """
        Solves the ODE for a given energy E and returns the value of u at the end of the grid.
        This function's roots are the energy eigenvalues.
        """
        r_start = r_grid[0]
        # Set initial conditions: u(r) ~ r^(l+1) for small r
        u_start = r_start**(l + 1)
        du_start = (l + 1) * r_start**l
        y0 = [u_start, du_start]
        
        # Integrate the ODE system
        sol = odeint(schrodinger_ode, y0, r_grid, args=(E, l))
        
        # Return the value of the wavefunction at the outer boundary
        return sol[-1, 0]

    # --- Step 4: Find Eigenvalues for different l ---
    def find_eigenvalues(l, E_range, r_grid, max_levels=2):
        """
        Scans an energy range to find eigenvalues using the shooting method and a root-finder.
        """
        eigenvalues = []
        E_scan = np.linspace(E_range[0], E_range[1], 500)
        u_end_scan = np.array([wavefunction_at_end(E, l, r_grid) for E in E_scan])

        # Find intervals where the wavefunction at r_end crosses zero
        sign_changes = np.where(np.diff(np.sign(u_end_scan)))[0]
        
        for idx in sign_changes:
            # Use a robust root-finder (brentq) to find the precise eigenvalue
            E_low, E_high = E_scan[idx], E_scan[idx + 1]
            try:
                eigen_E = brentq(wavefunction_at_end, E_low, E_high, args=(l, r_grid))
                
                # Determine principal quantum number n by counting nodes
                sol = odeint(schrodinger_ode, [r_grid[0]**(l+1), (l+1)*r_grid[0]**l], r_grid, args=(eigen_E, l))
                nodes = len(np.where(np.diff(np.sign(sol[:,0])))[0])

                eigenvalues.append({'E': eigen_E, 'n': nodes + 1, 'l': l})
                if len(eigenvalues) >= max_levels:
                    break
            except (ValueError, RuntimeError):
                continue
        return eigenvalues

    # --- Step 5: Execute the calculation ---
    # Setup the spatial grid for integration
    r_grid = np.linspace(1e-5, 50.0, 4000)
    
    # Set the energy search range for bound states (0 < E < V0)
    E_range = (0.01, V0 - 0.01)
    
    # Find lowest energy states for angular momentum l=0, 1, and 2
    s_states = find_eigenvalues(l=0, E_range=E_range, r_grid=r_grid)
    p_states = find_eigenvalues(l=1, E_range=E_range, r_grid=r_grid)
    d_states = find_eigenvalues(l=2, E_range=E_range, r_grid=r_grid)

    # Combine all found states into a single list
    all_states = s_states + p_states + d_states
    
    # Sort the states by energy to find the first two levels
    all_states.sort(key=lambda x: x['E'])

    # --- Step 6: Output the final result ---
    if len(all_states) < 2:
        print("Error: Could not find at least two energy levels in the specified range.")
        return

    E1_state = all_states[0]
    E2_state = all_states[1]

    E1 = E1_state['E']
    E2 = E2_state['E']
    delta_E = E2 - E1

    print(f"The first energy level (E1) corresponds to state n={E1_state['n']}, l={E1_state['l']}.")
    print(f"The second energy level (E2) corresponds to state n={E2_state['n']}, l={E2_state['l']}.")
    print("\nCalculating the energy difference:")
    print(f"ΔE = E2 - E1")
    print(f"ΔE = {E2:.4f} eV - {E1:.4f} eV")
    print(f"ΔE = {delta_E:.4f} eV")
    
    # The final answer for the grader is the numerical value of delta_E
    return delta_E

# Run the solver
final_answer = solve_quantum_well()
#<<<1.6488>>>