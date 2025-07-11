import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_critical_potential():
    """
    Calculates the critical chemical potential for holographic condensation
    in Einstein-Gauss-Bonnet gravity using a numerical shooting method.
    """
    # Parameters from the problem description
    lambda_gb = 0.1
    m2 = -3.0  # Corresponds to a scalar operator of dimension Delta=3

    # We work at a very low temperature, which corresponds to a large z_H (horizon radius).
    # The result for mu_c is independent of z_H for large z_H.
    z_H = 50.0

    # The metric function f(u) and its derivative, where u = z / z_H
    def f(u, lambda_gb):
        # Avoid u=1 singularity and sqrt of negative numbers
        if u >= 1.0:
            return 0.0
        # Check for causality constraint on lambda_gb
        if 1 - 4 * lambda_gb < 0:
             raise ValueError("Gauss-Bonnet coupling lambda_gb violates causality.")
        return (1 - np.sqrt(1 - 4 * lambda_gb * (1 - u**4))) / (2 * lambda_gb)

    def f_prime(u, lambda_gb):
        if u >= 1.0:
            # Analytical limit at u=1
            return 4.0
        if u == 0:
            return 0.0
        # Check for causality constraint on lambda_gb
        if 1 - 4 * lambda_gb < 0:
             raise ValueError("Gauss-Bonnet coupling lambda_gb violates causality.")
        return (4 * u**3) / np.sqrt(1 - 4 * lambda_gb * (1 - u**4))

    # The ODE for the scalar field perturbation psi at the critical point
    # We solve for y = [psi, psi_prime]
    def ode_system(u, y, eigenvalue, params):
        psi, psi_prime = y
        lambda_gb_val, m2_val = params
        
        f_u = f(u, lambda_gb_val)
        fp_u = f_prime(u, lambda_gb_val)

        if u == 0 or f_u == 0:
            return [psi_prime, 0]

        # This is the linearized equation for psi(u)
        psi_double_prime = -(fp_u / f_u - 3 / u) * psi_prime - \
                           ( (eigenvalue**2 * (1 - u**2)**2) / f_u**2 - m2_val / (u**2 * f_u) ) * psi
        
        return [psi_prime, psi_double_prime]

    # Asymptotic analysis near the boundary (u -> 0) gives the expected scaling exponent
    N_gb_sq = (1 - np.sqrt(1 - 4 * lambda_gb)) / (2 * lambda_gb)
    delta_plus = 2 + np.sqrt(4 + m2 / N_gb_sq)

    # Objective function for the root finder.
    # It returns the error from the desired boundary condition at u=0.
    def objective_function(eigenvalue):
        # Initial conditions determined by regularity at the horizon (u=1)
        u_start = 1.0 - 1e-5
        psi_start = 1.0  # Arbitrary normalization
        psi_prime_start = (m2 / f_prime(1.0, lambda_gb)) * psi_start
        
        y_start = [psi_start, psi_prime_start]
        
        # Integration range
        u_span = [u_start, 1e-5]

        # Solve the ODE
        sol = solve_ivp(
            ode_system, 
            u_span, 
            y_start, 
            args=([eigenvalue], [lambda_gb, m2]), 
            dense_output=True,
            method='RK45'
        )
        
        # Evaluate the solution at the final integration point (close to boundary)
        u_final = u_span[1]
        psi_final, psi_prime_final = sol.sol(u_final)
        
        # Check if the solution satisfies the boundary condition psi ~ u^delta_plus
        # This is equivalent to u * psi'/psi -> delta_plus for u->0
        error = (u_final * psi_prime_final / psi_final) - delta_plus
        return error

    # We are solving for the eigenvalue Lambda = mu_c * z_H
    # Use a root finding algorithm to find the value of Lambda that makes the objective function zero.
    # We need to provide a bracket [a, b] where objective_function(a) and objective_function(b) have opposite signs.
    # The value is known to be in the hundreds for z_H=50.
    try:
        sol = root_scalar(objective_function, bracket=[150, 250], method='bisect')
        Lambda = sol.root
        mu_c = Lambda / z_H
        print("The equation for the critical chemical potential, mu_c, has been solved numerically.")
        print(f"For a Gauss-Bonnet coupling lambda_GB = {lambda_gb}, the critical chemical potential is:")
        # We present the final value as a simple equation as requested.
        print(f"mu_c = {mu_c:.4f}")
    except ValueError as e:
        print(f"Could not find a root in the given bracket. The error was: {e}")
        print("This might happen if the initial bracket is incorrect or the problem is ill-defined for the given parameters.")

if __name__ == '__main__':
    solve_critical_potential()