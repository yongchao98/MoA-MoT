import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

def solve_quantum_well_problem():
    """
    Calculates the energy difference between the first and second energy levels
    of a particle in a given 3D potential well using the shooting method.
    """
    # Suppress warnings from the solver for cleaner output
    warnings.filterwarnings("ignore", category=UserWarning)

    # --- Constants and Parameters ---
    V0 = 15.0  # Potential energy at center in eV
    R = 3.0    # Well radius in nm
    M_PARTICLE = 9.11e-31 # kg
    
    # Use physical constants in eV and nm units for convenience
    HBARC = 197.327 # Planck's constant * speed of light (eV nm)
    C_LIGHT = 299792458 # m/s
    M_EV = M_PARTICLE * C_LIGHT**2 / 1.60218e-19 # Mass in eV/c^2
    
    # The crucial factor in the Schrodinger equation, in eV nm^2
    HBAR2_2M = (HBARC**2) / (2 * M_EV)

    # --- Potential Energy Function V(r) ---
    def potential(r, l):
        """
        Calculates the effective potential V_eff(r) = V(r) + V_centrifugal(r)
        for a given radius r and angular momentum quantum number l.
        """
        if r <= 0:
            return float('inf')

        # Calculate V(r) from the given V^2(r)
        if r < R:
            # For r < R, V^2(r) = V0 + W(exp(r - R))
            # Use .real to ensure the result is a real number
            V_r_sq = V0 + lambertw(np.exp(r - R)).real
        else: # r >= R
            # For r >= R, V^2(r) = V0 * (1 - (r/R)^-2)
            arg = 1 - (R / r)**2
            # Prevent numerical errors leading to sqrt of a negative number
            V_r_sq = V0 * max(0, arg)
            
        V_r = np.sqrt(V_r_sq)
        
        # Add the centrifugal term for l > 0
        V_centrifugal = HBAR2_2M * l * (l + 1) / r**2
        return V_r + V_centrifugal

    # --- ODE System for the Radial Schrödinger Equation ---
    def schrodinger_ode(r, y, E, l):
        """
        Represents the Schrodinger equation as a system of first-order ODEs.
        y = [u, u'] where u is the radial wavefunction.
        """
        u, u_prime = y
        du_dr = u_prime
        d_u_prime_dr = -(E - potential(r, l)) / HBAR2_2M * u
        return [du_dr, d_u_prime_dr]

    # --- Shooting Function ---
    def shoot(E, l, r_start, r_max):
        """
        Integrates the Schrodinger ODE for a given energy E and returns the
        wavefunction's value at the outer boundary r_max.
        """
        # Set initial conditions based on the behavior of u(r) near r=0
        if l == 0:
            y0 = [r_start, 1.0]  # For l=0, u(r) ~ r
        else:
            y0 = [r_start**(l + 1), (l + 1) * r_start**l] # For l>0, u(r) ~ r^(l+1)
        
        # Numerically solve the ODE
        sol = solve_ivp(
            schrodinger_ode,
            [r_start, r_max],
            y0,
            args=(E, l),
            t_eval=[r_max]
        )
        
        if sol.status != 0:
            return np.inf  # Return infinity if integration fails
            
        return sol.y[0, -1]

    # --- Eigenvalue Finder ---
    def find_eigenvalues(l, num_eigenvalues):
        """
        Finds the first 'num_eigenvalues' for a given l by finding the roots
        of the shooting function.
        """
        # Integration boundaries
        r_start = 1e-6 # A small number to avoid singularity at r=0
        r_max = 25.0   # A large radius where the wavefunction should decay

        # Energy search range for bound states [V_min, V_max_at_infinity]
        E_min = 0.01
        E_max = np.sqrt(V0)

        eigenvalues = []
        
        # Scan the energy range to find brackets where the solution changes sign
        energies = np.linspace(E_min, E_max - 0.01, 1000)
        u_final = np.array([shoot(E, l, r_start, r_max) for E in energies])
        
        sign_changes = np.where(np.diff(np.sign(u_final)))[0]
        
        for idx in sign_changes:
            if len(eigenvalues) >= num_eigenvalues:
                break
            
            E1_bracket, E2_bracket = energies[idx], energies[idx + 1]
            try:
                # Use a root-finder to get the precise eigenvalue
                eigenvalue = brentq(shoot, E1_bracket, E2_bracket, args=(l, r_start, r_max))
                eigenvalues.append(eigenvalue)
            except (ValueError, RuntimeError):
                # Ignore brackets where the root-finder fails
                continue
                
        return eigenvalues

    # --- Main Calculation ---
    
    # Find the ground state energy E1 (1s state)
    E_1s_levels = find_eigenvalues(l=0, num_eigenvalues=1)
    if not E_1s_levels:
        print("Error: Ground state (1s) not found. Cannot proceed.")
        return
    E1 = E_1s_levels[0]

    # Find the energy of the 2s state (first excited l=0 state)
    E_s_levels = find_eigenvalues(l=0, num_eigenvalues=2)
    E_2s = E_s_levels[1] if len(E_s_levels) > 1 else float('inf')
    
    # Find the energy of the 1p state (lowest l=1 state)
    E_1p_levels = find_eigenvalues(l=1, num_eigenvalues=1)
    E_1p = E_1p_levels[0] if E_1p_levels else float('inf')

    # The second energy level E2 is the minimum of the next available states
    E2 = min(E_2s, E_1p)
    if E2 == float('inf'):
        print("Error: Second energy level could not be determined.")
        return

    # Calculate the energy difference
    delta_E = E2 - E1

    # --- Output the final result ---
    print(f"The first energy level E1 (1s state) is: {E1:.6f} eV")
    print(f"The energy of the 2s state is: {E_2s:.6f} eV")
    print(f"The energy of the 1p state is: {E_1p:.6f} eV")
    print(f"The second energy level E2 is the lower of these two: {E2:.6f} eV")
    print("\nFinal calculation:")
    # The final equation as requested
    print(f"ΔE = E2 - E1 = {E2:.6f} - {E1:.6f} = {delta_E:.6f} eV")
    
    # Print the final answer in the required format
    print(f"\n<<<{delta_E:.6f}>>>")

# Execute the main function
solve_quantum_well_problem()