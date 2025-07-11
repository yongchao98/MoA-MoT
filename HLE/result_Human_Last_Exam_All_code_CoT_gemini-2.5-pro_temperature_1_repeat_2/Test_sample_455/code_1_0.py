import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

def solve_quantum_well_problem():
    """
    This function calculates the energy difference between the first and second
    energy levels of a quantum particle in the potential well defined in the problem.
    It uses a numerical shooting method to solve the radial Schrödinger equation.
    """
    
    # Step 1: Define physical constants and numerical parameters
    V0 = 15.0  # eV
    R = 3.0    # nm
    # h-bar^2 / (2 * m_e) in units of eV * nm^2
    HBAR_SQ_OVER_2M = 0.0380998  # eV*nm^2

    # Numerical parameters for the solver
    R_START = 1e-6    # A small radius to start integration, avoiding the r=0 singularity
    R_MAX = 20 * R    # A sufficiently large radius to check the boundary condition at infinity
    E_MIN_SEARCH = 0.01  # Search for energies above the potential minimum (0 eV)
    E_MAX_SEARCH = V0 - 0.01 # Search for bound states below the asymptotic potential (V0 = 15 eV)

    # Step 2: Define the potential energy function U(r) based on the given V^2(r)
    def potential_U(r):
        """Calculates the potential energy U(r) in eV."""
        if r < R:
            # For r < R, U(r) = V0 + W(exp(r - R)).
            # np.real is used to ensure the output is real, as lambertw can return complex numbers
            # for negative arguments (though the argument here is always positive).
            return V0 + np.real(lambertw(np.exp(r - R)))
        else:
            # For r >= R, U(r) = V0 * (1 - (r/R)^-2)
            if r == 0: return np.inf # Should not be called in this branch, but as a safeguard.
            return V0 * (1.0 - (R / r)**2)

    # Step 3: Define the effective potential, which includes the centrifugal barrier
    def U_eff(r, l):
        """Calculates the effective potential U_eff(r, l) for a given angular momentum l."""
        if r == 0:
            return np.inf
        centrifugal_term = l * (l + 1) * HBAR_SQ_OVER_2M / r**2
        return potential_U(r) + centrifugal_term

    # Step 4: Define the ODE system for the radial Schrödinger equation
    def radial_ode(r, y, E, l):
        """
        Defines the system of ODEs for u(r). y is a vector [u, u'].
        Returns dy/dr = [u', u''].
        """
        u, du_dr = y
        if r == 0:
            return [du_dr, 0]
        # This is the rearranged Schrödinger equation: u'' = (2m/hbar^2) * (U_eff - E) * u
        d2u_dr2 = (U_eff(r, l) - E) / HBAR_SQ_OVER_2M * u
        return [du_dr, d2u_dr2]

    # Step 5: Implement the core of the shooting method
    def shoot_and_get_final_u(E, l):
        """
        Solves the ODE for a given energy E and returns the wavefunction value at R_MAX.
        It handles the potential discontinuity at r=R by splitting the integration.
        """
        # Suppress benign warnings from the integrator for cleaner output
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Initial conditions near r=0, where u(r) is proportional to r^(l+1)
            u0 = [R_START**(l + 1), (l + 1) * R_START**l]

            # Integrate from R_START to the discontinuity point at R
            sol1 = solve_ivp(
                fun=radial_ode, t_span=[R_START, R], y0=u0, args=(E, l),
                method='RK45', atol=1e-9, rtol=1e-9
            )
            if sol1.status != 0: return np.inf # Return large value if integration fails

            # Use the wavefunction value at R as the initial condition for the next segment
            u_at_R = sol1.y[:, -1]

            # Integrate from R to R_MAX
            sol2 = solve_ivp(
                fun=radial_ode, t_span=[R, R_MAX], y0=u_at_R, args=(E, l),
                method='RK45', atol=1e-9, rtol=1e-9
            )
            if sol2.status != 0: return np.inf
            
            # For an eigenvalue, the wavefunction should be zero at infinity.
            # We return its value at R_MAX; the root-finder will find when this is zero.
            return sol2.y[0, -1]

    # Step 6: Create a function to find the n-th eigenvalue for a given l
    def find_eigenvalue(n_prime, l):
        """
        Finds the n'-th energy eigenvalue for a given angular momentum l.
        n_prime = 1 for the lowest energy state of that l, 2 for the next, etc.
        """
        # Scan over the allowed energy range to find intervals where eigenvalues might exist.
        E_vals = np.linspace(E_MIN_SEARCH, E_MAX_SEARCH, 300)
        f_vals = np.array([shoot_and_get_final_u(E, l) for E in E_vals])
        
        # Find indices where u(R_MAX) changes sign, indicating a root.
        sign_changes = np.where(np.diff(np.sign(f_vals)))[0]
        
        # Check if we found enough roots for the requested state
        if len(sign_changes) < n_prime:
            print(f"Error: Could not find state n'={n_prime} for l={l}. Try adjusting search parameters.")
            return None
        
        # Select the interval for the desired root (the n_prime-th sign change)
        interval_index = sign_changes[n_prime - 1]
        E_low, E_high = E_vals[interval_index], E_vals[interval_index + 1]
        
        # Use a precise root-finder (Brent's method) to find the eigenvalue in this interval.
        try:
            eigen_E = brentq(lambda E: shoot_and_get_final_u(E, l), E_low, E_high, xtol=1e-12)
            return eigen_E
        except ValueError:
            print(f"Root finding failed in interval [{E_low}, {E_high}].")
            return None

    # Step 7: Find the specific energy levels required by the problem
    # E1 is the ground state: l=0, n'=1 (lowest energy)
    E1 = find_eigenvalue(n_prime=1, l=0)

    # Candidates for the second energy level E2 are the next two lowest-energy states
    # First candidate: the first excited s-state (l=0, n'=2)
    E_2s = find_eigenvalue(n_prime=2, l=0)
    
    # Second candidate: the lowest p-state (l=1, n'=1)
    E_1p = find_eigenvalue(n_prime=1, l=1)
    
    if E1 is None or E_2s is None or E_1p is None:
        print("Calculation failed: Not all required energy levels could be found.")
        return

    # The second energy level, E2, is the lower of the two candidates
    E2 = min(E_2s, E_1p)

    # Step 8: Calculate the energy difference and print the final result
    delta_E = E2 - E1

    print(f"The first energy level (ground state) is E1 = {E1:.5f} eV.")
    print(f"The second energy level is E2 = {E2:.5f} eV.")
    print("The energy difference is Delta_E = E2 - E1, which is calculated as:")
    print(f"{E2:.5f} eV - {E1:.5f} eV = {delta_E:.5f} eV")
    
    print(f"\n<<<"+f"{delta_E:.5f}>>>")

# Execute the main function to solve the problem
solve_quantum_well_problem()