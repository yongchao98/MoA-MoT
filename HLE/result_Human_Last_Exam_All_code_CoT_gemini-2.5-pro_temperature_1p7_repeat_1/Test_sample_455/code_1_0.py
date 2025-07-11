import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

def solve_quantum_well_problem():
    """
    This script calculates the energy difference between the first and second
    energy levels of a particle in a specified potential well.
    """
    # --- 1. Constants and Parameters (in eV and nm) ---
    V0 = 15.0  # eV
    R = 3.0    # nm
    # hbar^2 / (2 * mass) in eV*nm^2 for an electron
    HBAR2_2M = (197.327**2) / (2 * 511000)

    # The potential approaches sqrt(V0) at infinity. Bound states exist for E < sqrt(V0).
    V_inf = np.sqrt(V0)

    # --- 2. Potential Energy Function ---
    def V(r):
        """Potential energy V(r) in eV."""
        # Ensure r is not exactly zero to avoid math errors
        if r <= 0:
            return np.inf
        if r < R:
            # For r < R, V(r) = sqrt(V0 + W(exp(r-R)))
            # The lambertw function returns a complex number, but for a positive real
            # argument, the principal branch (k=0) is real.
            return np.sqrt(V0 + lambertw(np.exp(r - R)).real)
        else:
            # For r >= R, V(r) = sqrt(V0 * (1 - (R/r)^2))
            # The argument to sqrt is always non-negative for r>=R
            return np.sqrt(V0 * (1 - (R/r)**2))

    # --- 3. ODE Solver (Shooting Method) ---
    def get_wavefunction_at_rmax(E, l, r_max):
        """
        Solves the radial Schrodinger equation for a given energy E and angular momentum l.
        Returns the value of the radial wavefunction u(r) at a large radius r_max.
        """
        # A valid bound state must have energy above the potential minimum (0)
        # and below the asymptotic potential at infinity.
        if not (0 < E < V_inf):
            return np.inf

        def V_eff(r):
            # Effective potential includes the centrifugal barrier
            centrifugal_term = HBAR2_2M * l * (l + 1) / (r**2)
            return V(r) + centrifugal_term

        # The ODE: u''(r) = (V_eff(r) - E) / (hbar^2/2m) * u(r)
        def schrodinger_ode(r, y):
            u, u_prime = y
            r_safe = r if r > 1e-9 else 1e-9  # Avoid division by zero
            d2u_dr2 = (V_eff(r_safe) - E) / HBAR2_2M * u
            return [u_prime, d2u_dr2]

        # Integrate from r_min to r_max.
        r_min = 1e-5

        # Set initial conditions based on the behavior of u(r) near the origin.
        if l == 0:
            y0 = [r_min, 1.0]  # u(r) ~ r for l=0
        else:
            y0 = [r_min**(l + 1), (l + 1) * r_min**l]  # u(r) ~ r^(l+1) for l>0

        # Numerically integrate the ODE
        sol = solve_ivp(schrodinger_ode, [r_min, r_max], y0, dense_output=True, atol=1e-9, rtol=1e-9)
        u_at_rmax = sol.sol(r_max)[0]
        return u_at_rmax

    # --- 4. Eigenvalue Finder ---
    def find_eigenvalues(l, E_min, E_max, r_max, num_levels_to_find):
        """Finds the first `num_levels_to_find` eigenvalues for a given l."""
        eigenvalues = []
        # Scan for energies where the wavefunction at r_max changes sign
        energy_steps = np.linspace(E_min, E_max, 700)
        u_prev = get_wavefunction_at_rmax(energy_steps[0], l, r_max)

        for i in range(1, len(energy_steps)):
            E_current = energy_steps[i]
            u_current = get_wavefunction_at_rmax(E_current, l, r_max)
            # A sign change indicates a root (eigenvalue) is bracketed
            if np.sign(u_current) != np.sign(u_prev):
                E_low, E_high = energy_steps[i - 1], energy_steps[i]
                try:
                    # Find the root precisely using Brent's method
                    eigen_E = brentq(lambda E: get_wavefunction_at_rmax(E, l, r_max), E_low, E_high)
                    eigenvalues.append(eigen_E)
                    if len(eigenvalues) == num_levels_to_find:
                        break
                except (ValueError, RuntimeError):
                    pass # Brentq can fail in some edge cases
            u_prev = u_current
        return eigenvalues

    # --- 5. Main Calculation ---
    # The integration radius r_max must be large enough for the wavefunction to decay.
    r_max = 50.0  # nm

    # Energy search range must be within the potential well: (0, V_inf)
    E_search_min = 0.01
    E_search_max = V_inf - 0.01

    # Find the lowest energy states for l=0, 1, and 2 to determine E1 and E2.
    energies_l0 = find_eigenvalues(l=0, E_min=E_search_min, E_max=E_search_max, r_max=r_max, num_levels_to_find=2)
    energies_l1 = find_eigenvalues(l=1, E_min=E_search_min, E_max=E_search_max, r_max=r_max, num_levels_to_find=1)
    energies_l2 = find_eigenvalues(l=2, E_min=E_search_min, E_max=E_search_max, r_max=r_max, num_levels_to_find=1)
    
    # Collect all found energies and sort them to find the true E1 and E2
    all_energies = sorted(energies_l0 + energies_l1 + energies_l2)

    if len(all_energies) < 2:
        print("Could not find at least two energy levels. Try adjusting r_max or search range.")
        return

    E1 = all_energies[0]
    E2 = all_energies[1]
    delta_E = E2 - E1

    # --- 6. Print the Final Result ---
    print(f"The ground state energy is E1 = {E1:.4f} eV.")
    print(f"The first excited state energy is E2 = {E2:.4f} eV.")
    print("\n--- Final Calculation ---")
    print(f"The energy difference is Delta E = E2 - E1")
    print(f"Delta E = {E2:.4f} eV - {E1:.4f} eV = {delta_E:.4f} eV")

    # Output the final answer in the requested format
    print(f"\n<<<{delta_E:.4f}>>>")

if __name__ == '__main__':
    # Suppress warnings that can occur during numerical integration of stiff ODEs
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    solve_quantum_well_problem()