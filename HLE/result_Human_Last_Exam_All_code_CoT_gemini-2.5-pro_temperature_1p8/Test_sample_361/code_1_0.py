import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a
    holographic D3/D7 model with an Einstein-Gauss-Bonnet background.
    """
    # The Gauss-Bonnet coupling parameter.
    lambda_gb = 0.1

    def f_metric(z, lgb):
        """
        Defines the Einstein-Gauss-Bonnet blackening function f(z).
        z=1 is the horizon, z=0 is the boundary. We set the horizon radius r_h=1.
        """
        # Ensure the argument of the square root is non-negative for numerical stability.
        sqrt_arg = 1 - 4 * lgb * (1 - z**4)
        if isinstance(sqrt_arg, np.ndarray):
            sqrt_arg[sqrt_arg < 0] = 0
        elif sqrt_arg < 0:
            sqrt_arg = 0
        
        # Note: For lgb -> 0, this expression correctly limits to 1 - z^4.
        if abs(lgb) < 1e-9:
            return 1 - z**4
            
        return (1 - np.sqrt(sqrt_arg)) / (2 * lgb)

    def ode_system(z, Y, mu, lgb):
        """
        Defines the system of first-order ODEs for the scalar field fluctuation psi.
        Y = [psi, psi_prime]
        """
        psi, psi_prime = Y
        
        # To avoid division by zero right at the horizon z=1 where f(z)=0,
        # we note that the full expression for psi'' is regular.
        # The solver will typically not step on z=1 exactly, but this is a safeguard.
        if z == 1.0:
            return [psi_prime, 0.0]
            
        f_val = f_metric(z, lgb)
        
        # This expression for f'/f is more stable numerically than direct computation.
        f_prime_over_f = (-4 * z**3) / (f_val * (1 - 2 * lgb * f_val))

        # The second-order ODE for psi(z)
        psi_double_prime = -(f_prime_over_f - 1/z) * psi_prime - (4 * mu**2 / f_val**2) * psi
        
        return [psi_prime, psi_double_prime]

    def objective_function(mu):
        """
        Objective function for the shooting method.
        It returns the value of psi at the boundary for a given mu.
        The root of this function gives the critical chemical potential mu_c.
        """
        # Integrate from just outside the horizon to near the boundary.
        z_start = 1.0 - 1e-5
        z_end = 1e-6

        # Set initial conditions for a regular solution at the horizon: psi(1)=1, psi'(1)=0.
        y0 = [1.0, 0.0]

        # Solve the initial value problem.
        sol = solve_ivp(
            ode_system,
            (z_start, z_end),
            y0,
            args=(mu, lambda_gb),
            method='RK45'
        )
        
        # The value of psi at the end of the integration interval.
        psi_at_boundary = sol.y[0, -1]
        
        return psi_at_boundary

    # Find the root of the objective function to get mu_c.
    # The range for mu is estimated based on known results in similar models.
    try:
        mu_c = brentq(objective_function, 0.8, 1.2, xtol=1e-6)
        print(f"The equation for the scalar fluctuation profile psi(z) has been solved numerically.")
        print(f"The parameters for the equation are:")
        print(f"Gauss-Bonnet coupling (lambda_GB): {lambda_gb}")
        print(f"\nThe value of the critical chemical potential (mu_c) is:")
        print(f"{mu_c:.8f}")

    except ValueError:
        print("Failed to find the root. The initial bracketing values might be incorrect.")

# Execute the calculation and print the result.
solve_critical_potential()