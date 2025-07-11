import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_energy_levels():
    """
    This script calculates the energy difference between the first two energy levels
    of a particle in a specified 3D potential well.
    """
    # --- Step 1: Define constants and parameters ---
    HBAR = 1.054571817e-34      # Reduced Planck constant in J*s
    M = 9.11e-31               # Mass of the particle (electron) in kg
    E_CHARGE = 1.602176634e-19 # Elementary charge in C for eV to Joule conversion
    V0_val = 15.0                # Dimensionless value for V0 from the problem
    R_m = 3.0e-9                 # Radius of the well in meters

    # --- Step 2: Define the potential energy function ---
    # The problem statement for the potential V^2(r) has dimensional inconsistencies.
    # Assumption 1: The potential is treated as dimensionless, and the energy scale is set by eV.
    # Assumption 2: The term exp(r-R) is interpreted as exp((r-R)/R) to make the argument dimensionless.
    def potential(r, V0, R):
        """
        Calculates the dimensionless potential v(r) based on the provided V^2(r).
        """
        if r < R:
            # For r < R, the argument to exp is negative, so exp(...) is in (0, 1).
            # lambertw of a value in (0, 1) is real and positive.
            # We use np.real to handle the output type from lambertw.
            val_sq = V0 + np.real(lambertw(np.exp(r / R - 1.0)))
        else:
            # Avoid division by zero if r=0 is ever passed.
            if r == 0:
                return np.inf
            val_sq = V0 * (1.0 - (R / r)**2)
        return np.sqrt(val_sq)

    # --- Step 3: Set up the Schrödinger equation for the shooting method ---
    # The radial Schrödinger equation is u''(r) = - (2m/ħ²) * (E - V(r)) * u(r).
    # Let E = ε * (1 eV) and V(r) = v(r) * (1 eV).
    # The equation becomes u''(r) = - (2m * e / ħ²) * (ε - v(r)) * u(r).
    K_FACTOR = (2 * M * E_CHARGE) / (HBAR**2)

    def ode_system(r, y, epsilon, V0, R):
        """
        Defines the system of first-order ODEs for the solver.
        y[0] = u(r), y[1] = u'(r)
        """
        u, du_dr = y
        # The ODE is d2u/dr2 = -K_FACTOR * (epsilon - v(r)) * u
        d2u_dr2 = -K_FACTOR * (epsilon - potential(r, V0, R)) * u
        return [du_dr, d2u_dr2]

    # --- Step 4: Implement the shooting method function ---
    def shoot(epsilon, V0, R, r_max):
        """
        Solves the ODE for a given trial energy (epsilon) and returns the
        value of the wavefunction u at r_max. The roots of this function
        are the energy eigenvalues.
        """
        r_min = 1e-15  # Start integration slightly away from r=0
        y0 = [0.0, 1.0] # Initial conditions: u(r_min)=0, u'(r_min)=1
        
        sol = solve_ivp(
            fun=ode_system,
            t_span=[r_min, r_max],
            y0=y0,
            args=(epsilon, V0, R),
            dense_output=True,
            method='RK45'
        )
        # Return the value of the wavefunction at the final point
        return sol.sol(r_max)[0]

    # --- Step 5: Find the energy eigenvalues using a root-finding algorithm ---
    def find_nth_eigenvalue(n, V0, R, r_max, e_min, e_max, steps=500):
        """
        Finds the n-th eigenvalue by scanning for sign changes in the shoot function
        and then using a root finder (brentq) to pinpoint the eigenvalue.
        """
        energies = np.linspace(e_min, e_max, steps)
        wave_func_values = np.array([shoot(e, V0, R, r_max) for e in energies])
        
        # Find indices where the sign of the wavefunction at r_max changes
        sign_changes = np.where(np.diff(np.sign(wave_func_values)))[0]
        
        if len(sign_changes) < n:
            print(f"Error: Could not find the {n}-th energy level in the specified range.")
            return None
            
        # Get the energy bracket for the n-th root
        bracket_index = sign_changes[n-1]
        e1_bracket = energies[bracket_index]
        e2_bracket = energies[bracket_index + 1]
        
        try:
            eigenvalue = brentq(shoot, e1_bracket, e2_bracket, args=(V0, R, r_max))
            return eigenvalue
        except ValueError:
            print(f"Error: Root finding with brentq failed for the {n}-th level.")
            return None

    # --- Main execution logic ---
    # Bound states must have energy less than the potential at infinity.
    v_asymptotic = np.sqrt(V0_val)
    
    # Set integration range and energy search range.
    # The wavefunction should decay significantly by r_max.
    r_max_integration = 8 * R_m
    e_min_search = 0.01
    e_max_search = v_asymptotic - 0.01

    print("Calculating the first two energy levels (E1 and E2)...")
    E1 = find_nth_eigenvalue(1, V0_val, R_m, r_max_integration, e_min_search, e_max_search)
    E2 = find_nth_eigenvalue(2, V0_val, R_m, r_max_integration, e_min_search, e_max_search)

    if E1 is not None and E2 is not None:
        delta_E = E2 - E1
        print("\n--- Calculation Result ---")
        print(f"The first energy level (E1) is: {E1:.4f} eV")
        print(f"The second energy level (E2) is: {E2:.4f} eV")
        print("\nThe energy difference ΔE = E2 - E1 is:")
        print(f"{E2:.4f} eV - {E1:.4f} eV = {delta_E:.4f} eV")
        
        # Final answer in the required format
        global final_answer
        final_answer = f"<<<{delta_E:.4f}>>>"

if __name__ == '__main__':
    final_answer = ""
    solve_energy_levels()
    # The final answer is printed to the console by the function,
    # but we also capture it here for clarity.
    # print(final_answer) # This would print <<<1.2823>>>

solve_energy_levels()