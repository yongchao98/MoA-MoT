import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a 5D
    Einstein-Gauss-Bonnet holographic model.
    """
    # --- 1. Define Model Parameters ---
    LAMBDA_GB = 0.1  # Gauss-Bonnet coupling
    M2 = -4.0        # Scalar mass squared for operator dimension Delta=2
    Q = 1.0          # Scalar charge
    ZH = 1.0         # Set horizon radius to 1 (this sets the temperature scale)

    # --- 2. Define Background Functions ---
    def f(z, lambda_gb):
        """EGB black hole metric function f(z)."""
        discriminant = 1 - 4 * lambda_gb * (1 - z**4)
        return (1.0 - np.sqrt(discriminant)) / (2.0 * lambda_gb)

    def f_prime(z, lambda_gb):
        """Derivative of f(z) with respect to z."""
        discriminant = 1 - 4 * lambda_gb * (1 - z**4)
        return 4 * z**3 / np.sqrt(discriminant)

    def phi(z, mu):
        """Gauge field profile phi(z) = A_t(z)."""
        return mu * (1 - z**2)

    # --- 3. Define the ODE System ---
    def ode_system(z, y, mu, lambda_gb):
        """
        Defines the second-order ODE for the scalar field psi as a system
        of first-order ODEs. Y = [psi, psi_prime].
        """
        psi, psi_prime = y
        
        f_val = f(z, lambda_gb)
        # Handle numerical precision issues very close to the horizon
        if np.isclose(f_val, 0):
            return [psi_prime, 1e18]

        f_p_val = f_prime(z, lambda_gb)
        phi_val = phi(z, mu)

        # The ODE is: psi'' + P(z)psi' + Q(z)psi = 0
        P_z = f_p_val / f_val - 3.0 / z
        Q_z = (Q**2 * phi_val**2) / (f_val**2) - M2 / (z**2 * f_val)
        
        psi_double_prime = -P_z * psi_prime - Q_z * psi
        return [psi_prime, psi_double_prime]

    # --- 4. Implement the Shooting Method ---
    def shoot(mu):
        """
        Shooting function that returns a residue. The root of this function
        is the critical chemical potential mu_c.
        """
        # Integrate from just outside the horizon to just outside the boundary
        epsilon_h = 1e-6  # Near horizon z = 1
        epsilon_b = 1e-5  # Near boundary z = 0
        z_span = [ZH - epsilon_h, epsilon_b]

        # Boundary condition at the horizon (z=1) from regularity
        # f'(1)psi'(1) - m^2 psi(1) = 0
        # With f'(1)=-4 and m^2=-4, this gives psi'(1) = psi(1).
        psi_init = 1.0
        psi_prime_init = 1.0
        y0 = [psi_init, psi_prime_init]

        # Solve the initial value problem
        sol = solve_ivp(
            ode_system, z_span, y0, args=(mu, LAMBDA_GB),
            method='RK45', dense_output=True, atol=1e-8, rtol=1e-8
        )
        
        # If solver fails, return a large value to guide root-finder
        if sol.status != 0:
            return 1e18

        # Extract solution at the boundary end of the integration
        z_end = epsilon_b
        psi_end, psi_prime_end = sol.sol(z_end)

        # Boundary condition at z=0 for condensation (zero source)
        # For m^2=-4 (Delta=2), the solution should behave as psi ~ z^2.
        # This implies that the quantity (2 * psi / z - psi') should be zero.
        residue = 2.0 * psi_end / z_end - psi_prime_end
        return residue

    # --- 5. Find the Root to get mu_c ---
    # We search for the root in a physically motivated bracket.
    bracket = [4.0, 8.0]
    try:
        result = root_scalar(shoot, bracket=bracket, method='brentq')
        if result.converged:
            mu_c = result.root
            print(f"Calculation for Gauss-Bonnet coupling lambda_GB = {LAMBDA_GB}")
            print(f"The critical chemical potential is:")
            print(mu_c)
        else:
            print("Root finding did not converge.")
    except ValueError:
        print(f"The shooting function did not have opposite signs at the bracket endpoints {bracket}.")
        print("Please try a different bracket.")

if __name__ == '__main__':
    solve_critical_potential()