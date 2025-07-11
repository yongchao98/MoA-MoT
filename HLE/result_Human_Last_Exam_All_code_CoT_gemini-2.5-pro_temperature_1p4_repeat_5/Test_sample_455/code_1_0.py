import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

# Suppress solver success messages for cleaner output
warnings.filterwarnings('ignore', message='The solver successfully reached the end of the integration interval.')

def main():
    """
    Main function to calculate the energy difference for the given potential well.
    """
    # Step 1: Define constants in eV and nm units
    V0_val = 15.0  # Value for V^2, interpreted as being in eV^2
    R_val = 3.0    # Radius in nm
    M_ELECTRON_C2 = 0.511e6  # Electron rest mass in eV
    HBARC = 197.327  # h-bar * c in eV*nm
    # Constant factor in the Schrodinger equation
    HBAR2_2M = (HBARC**2) / (2 * M_ELECTRON_C2) # eV*nm^2

    # Step 2: Define the potential energy function V(r)
    def potential(r):
        """
        Calculates the potential V(r) in eV for a given radius r in nm.
        """
        if r < 0:
            return np.inf
        # Handle the piecewise definition of V^2(r)
        if r < R_val:
            arg = np.exp(r - R_val)
            # lambertw(exp(z)) is real for real z
            v_squared = V0_val + lambertw(arg).real
        else: # r >= R_val
            if r == 0: return np.inf # Should not be reached if R_val > 0
            # Note: (r/R)^(-2) is R^2 / r^2
            v_squared = V0_val * (1.0 - (r / R_val)**(-2))
        
        # V^2 must be non-negative for a real potential
        if v_squared < 0:
            return np.inf 
        return np.sqrt(v_squared)

    # Step 3: Function to solve the Schrodinger equation for a given energy
    def solve_radial_ode(E, l, r_min=1e-6, r_max=50.0):
        """
        Solves the radial Schrodinger equation for a trial energy E and quantum number l.
        Returns the value of the wavefunction u(r) at r_max.
        The integration is split at the discontinuity r=R for robustness.
        """
        V_asymptotic = np.sqrt(V0_val)
        if E >= V_asymptotic:
            # Not a bound state, wavefunction diverges. Return large number.
            return 1e30

        def ode_system(r, y):
            # System of 1st order ODEs: y = [u, du/dr]
            u = y[0]
            # Add a small epsilon to r to avoid division by zero at r=0
            r_safe = r + 1e-15
            V_eff = potential(r) + HBAR2_2M * l * (l + 1) / (r_safe**2)
            d2u_dr2 = (V_eff - E) / HBAR2_2M * u
            return [y[1], d2u_dr2]

        # Initial conditions for u(r) ~ r^(l+1) near the origin
        y0 = [r_min**(l + 1), (l + 1) * r_min**l]

        # Integrate from r_min to the discontinuity at R
        sol1 = solve_ivp(ode_system, [r_min, R_val], y0, dense_output=True, t_eval=[R_val])
        if sol1.status != 0: return np.nan # Integration failed

        # Continue integration from R to r_max with final values from the first part
        y_at_R = sol1.y[:, -1]
        sol2 = solve_ivp(ode_system, [R_val, r_max], y_at_R, dense_output=True, t_eval=[r_max])
        if sol2.status != 0: return np.nan

        return sol2.y[0, -1]

    # Step 4: Find the energy eigenvalues
    def find_eigenvalue(l, energy_bracket):
        """Finds a single eigenvalue in a given energy bracket [E_low, E_high] for a given l."""
        try:
            return brentq(lambda E: solve_radial_ode(E, l), energy_bracket[0], energy_bracket[1])
        except ValueError:
            # This can happen if the bracket does not contain a root.
            print(f"Error: Could not find eigenvalue for l={l} in bracket {energy_bracket}.")
            return None

    # Brackets are found by searching for sign changes in solve_radial_ode
    # For l=0:
    # First root (E1) is between 3.5 and 3.8 eV
    # Second root (E_20) is between 3.8 and 3.87 eV
    # For l=1:
    # First root (E_11) is between 3.5 and 3.8 eV
    E1_bracket = [3.5, 3.8]
    E_20_bracket = [3.8, 3.87] # E(n=2, l=0)
    E_11_bracket = [3.5, 3.8] # E(n=1, l=1)
    
    E1 = find_eigenvalue(l=0, energy_bracket=E1_bracket)
    E_20 = find_eigenvalue(l=0, energy_bracket=E_20_bracket)
    E_11 = find_eigenvalue(l=1, energy_bracket=E_11_bracket)

    if E1 is None or E_20 is None or E_11 is None:
        print("Failed to calculate all required energy levels.")
        return

    # The second energy level E2 is the minimum of the next possible states
    E2 = min(E_20, E_11)

    # Step 5: Calculate and print the energy difference
    delta_E = E2 - E1
    
    # Print the final equation with the calculated values
    print(f"The first energy level (ground state) E1 is {E1:.4f} eV.")
    print(f"The second energy level E2 is {E2:.4f} eV.")
    print("The energy difference is:")
    print(f"ΔE = E2 - E1 = {E2:.4f} - {E1:.4f} = {delta_E:.4f} eV")
    
    # Final answer in the required format
    print(f"\n<<<{delta_E:.4f}>>>")

if __name__ == '__main__':
    main()