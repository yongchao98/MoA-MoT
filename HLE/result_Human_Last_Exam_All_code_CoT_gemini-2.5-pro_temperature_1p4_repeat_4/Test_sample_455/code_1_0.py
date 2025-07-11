import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import lambertw
from scipy.optimize import brentq

def solve_quantum_well_problem():
    """
    Solves the Schrodinger equation for the given potential to find the
    energy difference between the first two energy levels.
    """
    # --- 1. Constants and Parameters (in eV, nm) ---
    # Using units of eV for energy and nm for length simplifies the problem.
    m_electron_ev = 0.511e6  # Mass of electron in eV/c^2
    hbar_c = 197.327          # h-bar * c in eV*nm
    # The kinetic energy term constant hbar^2 / (2m)
    HBAR2_OVER_2M = hbar_c**2 / (2 * m_electron_ev)  # ~0.0381 eV*nm^2

    # Problem parameters
    V0_param_eV = 15.0  # The parameter V0 in the potential function, in eV
    R_nm = 3.0          # The radius of the well, in nm

    # --- 2. Potential Energy Function (in eV, nm) ---
    def potential_eV(r, V0, R):
        """Calculates potential V(r) in eV, with r and R in nm."""
        # This potential is discontinuous at r=R.
        # It forms a repulsive core for r<R and a well for r>=R.
        if r < R:
            # Assume r and R are their numerical values in nm for the exponential.
            arg = np.exp(r - R)
            v_sq = V0 + np.real(lambertw(arg))
        else:  # r >= R
            v_sq = V0 * (1.0 - (R / r)**2)
        
        # The potential squared should not be negative.
        if v_sq < 0:
            return np.inf  # Return a large number to signify a forbidden region.
        return np.sqrt(v_sq)

    # --- 3. Schrödinger Equation Solver (Shooting Method) ---
    def ode_system_ev(r, y, E, l):
        """Defines the system of ODEs for the radial Schrödinger equation."""
        u, du_dr = y
        
        # Centrifugal term, handle r=0 case
        if r == 0:
            V_centrifugal = 0
        else:
            V_centrifugal = HBAR2_OVER_2M * l * (l + 1) / (r**2)
            
        V_eff = potential_eV(r, V0_param_eV, R_nm) + V_centrifugal
        d2u_dr2 = (1 / HBAR2_OVER_2M) * (V_eff - E) * u
        return [du_dr, d2u_dr2]

    def shooting_function_ev(E, l, r_max):
        """
        Solves the ODE for a given energy E and returns u(r_max).
        The roots of this function are the energy eigenvalues.
        """
        r_min = 1e-6  # Start integration just above zero to avoid singularity.
        
        # Initial conditions for u(r) ~ r^(l+1) for small r
        u_min = r_min**(l + 1)
        du_dr_min = (l + 1) * r_min**l
        y_init = [u_min, du_dr_min]
        
        sol = solve_ivp(
            fun=lambda r, y: ode_system_ev(r, y, E, l),
            t_span=[r_min, r_max],
            y0=y_init,
            method='RK45', dense_output=True,
            atol=1e-10, rtol=1e-8
        )
        # Return the value of the wavefunction at the outer boundary.
        return sol.sol(r_max)[0]

    # --- 4. Eigenvalue Finding Function ---
    def find_eigenvalues_ev(l, num_states, E_min, E_max, r_max):
        """Searches for the first `num_states` eigenvalues for a given l."""
        eigenvalues = []
        # Create a fine grid of energies to scan for sign changes
        E_grid = np.linspace(E_min, E_max, 2000)
        
        f_prev = shooting_function_ev(E_grid[0], l, r_max)
        
        for i in range(1, len(E_grid)):
            E = E_grid[i]
            f_curr = shooting_function_ev(E, l, r_max)
            # A sign change indicates a root (eigenvalue) in the interval
            if np.sign(f_curr) != np.sign(f_prev):
                try:
                    # Use Brent's method for accurate root finding
                    eigen_E = brentq(
                        f=lambda e: shooting_function_ev(e, l, r_max),
                        a=E_grid[i-1], b=E
                    )
                    eigenvalues.append(eigen_E)
                    if len(eigenvalues) >= num_states:
                        break
                except ValueError:
                    # Brentq might fail if signs are the same due to precision, so we just continue
                    pass
            f_prev = f_curr
            
        return eigenvalues

    # --- 5. Main Calculation ---
    # The wavefunction should decay for r > R. We integrate far enough for it to be near zero.
    r_max_nm = 15.0  # 5 times the well radius R

    # Bound states must have energy E < V(infinity).
    V_infinity = np.sqrt(V0_param_eV)  # approx 3.873 eV
    E_search_min_eV = 0.01
    E_search_max_eV = V_infinity * 0.999 # Search just below the asymptote

    # Find the ground state energy E1 (1s state: l=0, lowest energy)
    E_1s_list = find_eigenvalues_ev(l=0, num_states=1, E_min=E_search_min_eV, E_max=E_search_max_eV, r_max=r_max_nm)
    if not E_1s_list:
        raise RuntimeError("Calculation failed: Ground state (1s) not found.")
    E1_eV = E_1s_list[0]

    # To find the first excited state E2, we compare the energies of the 2p and 2s states.
    # Find the lowest p-state energy (2p state: l=1, lowest energy for l=1)
    E_2p_list = find_eigenvalues_ev(l=1, num_states=1, E_min=E_search_min_eV, E_max=E_search_max_eV, r_max=r_max_nm)
    if not E_2p_list:
        raise RuntimeError("Calculation failed: 2p state not found.")
    E_2p_eV = E_2p_list[0]

    # Find the second s-state energy (2s state: l=0, second lowest energy for l=0)
    E_s_states = find_eigenvalues_ev(l=0, num_states=2, E_min=E_search_min_eV, E_max=E_search_max_eV, r_max=r_max_nm)
    if len(E_s_states) < 2:
        raise RuntimeError("Calculation failed: 2s state not found.")
    E_2s_eV = E_s_states[1]

    # The first excited state E2 is the minimum of the 2p and 2s energies.
    E2_eV = min(E_2p_eV, E_2s_eV)

    # Calculate the final energy difference
    delta_E_eV = E2_eV - E1_eV

    # --- 6. Output Results ---
    print("The final equation is: Delta E = E2 - E1")
    print(f"E1 (ground state energy): {E1_eV:.4f} eV")
    print(f"E2 (first excited state energy): {E2_eV:.4f} eV")
    print(f"Delta E (energy difference): {delta_E_eV:.4f} eV")
    
    return delta_E_eV

# Execute the calculation and store the final answer.
final_answer = solve_quantum_well_problem()
print(f"\n<<<Result value for verification>>>")
print(f"<<<{final_answer:.4f}>>>")
