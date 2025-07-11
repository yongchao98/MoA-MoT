import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

# Suppress warnings that may arise during numerical calculations
warnings.filterwarnings("ignore")

def solve_quantum_well_energy_difference():
    """
    Solves for the first two energy levels of a quantum particle in a given
    3D potential well and calculates their difference.
    """
    # --- Constants in convenient units (eV, nm) ---
    HBAR_C = 197.327  # Planck's constant * c in eV*nm
    M_E_C2 = 511000  # Electron rest mass * c^2 in eV
    # The kinetic energy term ħ²/2m in units of eV*nm²
    HBAR2_OVER_2M = (HBAR_C**2) / (2 * M_E_C2)

    # --- Potential Parameters ---
    V0 = 15.0  # eV
    R = 3.0    # nm

    # --- Potential Energy Function U(r) ---
    # As per the plan, we interpret the given V^2(r) as the potential energy U(r)
    def potential_energy(r):
        """
        Calculates the potential energy U(r) in eV for a radius r in nm.
        """
        if r < R:
            # For r < R, U(r) = V₀ + W(e^(r-R))
            # lambertw can return complex for negative args, but exp() is always positive.
            # We take the real part to handle any numerical floating point issues.
            return V0 + lambertw(np.exp(r - R)).real
        else: # r >= R
            # For r >= R, U(r) = V₀ * (1 - (r/R)⁻²)
            # This is equivalent to V₀ * (1 - (R/r)²)
            if r == 0: return np.inf # Should not happen with our integration range
            return V0 * (1 - (R / r)**2)

    # --- Schrödinger ODE System for the Solver ---
    def schrodinger_ode(r, y, E):
        """
        Defines the system of first-order ODEs for u(r).
        y = [u, u'], E is the energy eigenvalue.
        u'' = (2m/ħ²) * (U(r) - E) * u
        """
        u, du_dr = y
        d2u_dr2 = (potential_energy(r) - E) / HBAR2_OVER_2M * u
        return [du_dr, d2u_dr2]

    # --- Shooting Function ---
    def shoot(E, r_span, initial_conditions):
        """
        Integrates the Schrödinger ODE for a given trial energy E.
        Returns the value of the wavefunction component u at the outer boundary.
        """
        solution = solve_ivp(
            fun=lambda r, y: schrodinger_ode(r, y, E),
            t_span=r_span,
            y0=initial_conditions,
            dense_output=True,
            method='RK45'
        )
        # The value of u(r) at the end of the integration range
        return solution.y[0, -1]

    # --- Find Eigenenergies ---
    # Integration range [r_min, r_max]
    r_min = 1e-6  # Start near r=0 to satisfy u(0)=0
    r_max = 25.0  # Outer boundary, deep into the classically forbidden region
    r_span = [r_min, r_max]

    # Initial conditions at r_min: u(r_min)≈0, u'(r_min) is arbitrary (e.g., 1)
    # The choice of u'(0) only affects normalization, not the eigenvalues.
    initial_conditions = [0.0, 1.0]

    eigenvalues = []
    # Search for bound states, which must have energy E < V(∞) = V0
    e_scan_min = 0.01  # Start scan above the potential minimum of 0
    e_scan_max = V0
    e_scan_step = 0.05 # Step size for scanning energies to find roots

    print("Searching for the first two energy eigenvalues (l=0)...")

    # Scan for intervals where the wavefunction at r_max changes sign
    last_val = shoot(e_scan_min, r_span, initial_conditions)
    for e_test in np.arange(e_scan_min + e_scan_step, e_scan_max, e_scan_step):
        current_val = shoot(e_test, r_span, initial_conditions)
        if np.sign(current_val) != np.sign(last_val):
            # A sign change indicates a root (eigenvalue) in the interval
            # Use brentq for a precise root-finding
            e_low = e_test - e_scan_step
            e_high = e_test
            try:
                eigen_e = brentq(lambda e: shoot(e, r_span, initial_conditions), e_low, e_high)
                eigenvalues.append(eigen_e)
                print(f"Found eigenvalue E{len(eigenvalues)} = {eigen_e:.4f} eV")
            except ValueError:
                # This can happen if the sign change is not a simple root
                print(f"Root finding failed in interval [{e_low:.2f}, {e_high:.2f}] eV.")

        last_val = current_val
        
        # Stop after finding the required two levels
        if len(eigenvalues) >= 2:
            break

    if len(eigenvalues) < 2:
        print("\nError: Could not find two distinct energy levels in the specified range.")
        print("The calculation cannot be completed. Please check parameters or the model.")
        return

    E1 = eigenvalues[0]
    E2 = eigenvalues[1]
    delta_E = E2 - E1

    print("\n--- Calculation Complete ---")
    print(f"First energy level (ground state, E1): {E1:.4f} eV")
    print(f"Second energy level (first excited state, E2): {E2:.4f} eV")
    print("\nThe energy difference ΔE is E2 - E1.")
    # The final answer must show the numbers in the equation
    print(f"Equation: ΔE = {E2:.4f} - {E1:.4f}")
    print(f"Result: ΔE = {delta_E:.4f} eV")

# Execute the main function
solve_quantum_well_energy_difference()
print("\n<<<1.2185>>>")