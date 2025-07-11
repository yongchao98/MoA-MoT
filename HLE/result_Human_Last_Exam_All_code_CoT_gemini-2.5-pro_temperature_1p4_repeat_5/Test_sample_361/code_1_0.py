import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic s-wave superconductor
    in 5D Einstein-Gauss-Bonnet gravity.
    """
    # Suppress RuntimeWarning from integration near singularities
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    # Define model parameters
    LAMBDA_GB = 0.1  # Gauss-Bonnet coupling
    M2 = 0           # Mass squared of the scalar field (m^2=0 for a dim-4 operator)

    # --- Define the necessary functions for the ODE ---

    def f(z, lamb):
        """The metric function f(z) for EGB-AdS black hole."""
        # Use a safe z to avoid numerical errors at the boundaries
        z_safe = np.clip(z, 1e-12, 1.0 - 1e-12)
        sqrt_arg = 1 - 4 * lamb * (1 - z_safe**4)
        return (1 / (2 * lamb)) * (1 - np.sqrt(sqrt_arg))

    def df_dz(z, lamb):
        """The derivative of f(z) with respect to z."""
        z_safe = np.clip(z, 1e-12, 1.0 - 1e-12)
        sqrt_arg = 1 - 4 * lamb * (1 - z_safe**4)
        return (2 * z_safe**3) / np.sqrt(sqrt_arg)

    def odes(z, y, mu, lamb):
        """
        The system of ODEs for the scalar field psi.
        y[0] = psi, y[1] = psi'(z)
        """
        if z < 1e-9:  # At the boundary z=0, the equation is trivial
            return [y[1], 0]

        f_val = f(z, lamb)
        df_val = df_dz(z, lamb)

        if abs(f_val) < 1e-12:
             return [0, 0] # Return zero derivatives if f is too small

        fp_over_f = df_val / f_val
        
        # Coefficient of psi' in the standard form y'' + P(z)y' + Q(z)y = 0
        P_z = fp_over_f - 3 / z
        # Coefficient of psi
        Q_z = (mu**2 * (1 - z**2)**2) / (f_val**2) - M2 / (z**2 * f_val)
        
        dpsidz = y[1]
        d2psidz2 = -P_z * dpsidz - Q_z * y[0]
        
        return [dpsidz, d2psidz2]

    # --- Define the shooting function ---

    def shoot(mu, lamb=LAMBDA_GB):
        """
        Solves the ODE for a given mu and returns the value of psi at the boundary.
        The critical mu is found when psi(boundary) = 0.
        """
        # Integrate from just outside the horizon (z=1) to just before the boundary (z=0)
        z_start = 1 - 1e-6
        z_end = 1e-7

        # Regularity at the horizon implies initial conditions:
        # psi(1) = 1 (by normalization) and psi'(1) = 0.
        y0 = [1.0, 0.0]

        # Solve the ODE
        sol = solve_ivp(odes, [z_start, z_end], y0, args=(mu, lamb), dense_output=True, method='RK45')
        
        # Return the value of psi at the boundary
        psi_at_boundary = sol.sol(z_end)[0]
        return psi_at_boundary

    # --- Find the root and print the results ---

    print("To find the critical chemical potential, we solve the equations for a holographic")
    print("s-wave superconductor model in 5D Einstein-Gauss-Bonnet (EGB) gravity.")
    print("\nThe key parameters of the model are:")
    print(f"  - The Gauss-Bonnet coupling: lambda_GB = {LAMBDA_GB}")
    print(f"  - The scalar field mass squared: m^2 = {M2} (for a Delta=4 operator)")
    
    print("\nThe condensation is triggered when the chemical potential `mu` reaches a critical value `mu_c`.")
    print("We find `mu_c` by solving the linearized equation for the scalar field `psi`:")
    print("\n  psi''(z) + (f'(z)/f(z) - 3/z)psi'(z) + (mu^2 * (1 - z^2)^2 / f(z)^2)psi(z) = 0\n")
    print("where `z` is the radial coordinate (z=0 at boundary, z=1 at horizon), and f(z) is the EGB metric function.")

    print("We use a numerical shooting method to find the lowest `mu` that allows a solution where `psi(z=0) = 0`.")
    
    try:
        # Bracket the root. A quick check shows it's between 4 and 5.
        mu_c = brentq(shoot, 4.0, 5.0)
        print("\n--- Calculation Result ---")
        print(f"The critical chemical potential is calculated to be: mu_c = {mu_c:.4f}")

    except ValueError:
        print("\nCould not find a root in the given interval [4.0, 5.0]. The bracketing might be incorrect.")

    # Return the final value as per the requested format
    return f"<<<{mu_c:.4f}>>>"

# Run the calculation and print the final answer in the required format
final_answer = solve_critical_potential()
print(final_answer)