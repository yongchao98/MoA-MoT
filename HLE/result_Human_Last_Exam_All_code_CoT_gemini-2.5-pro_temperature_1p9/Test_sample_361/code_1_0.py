import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_for_mu_c():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using a numerical shooting method.

    The model parameters are fixed as follows:
    - Gauss-Bonnet coupling (lambda_gb): 0.1
    - Mass squared of the scalar field (m_squared): -3 (from a D3/D7 setup)
    - Temperature: 0 (extremal black hole)
    """

    # 1. Define the model parameters
    lambda_gb = 0.1
    m_squared = -3
    zh = 1.0  # Set horizon radius to 1. mu_c will be in units of zh.

    # 2. Define the EGB metric function f(z) and its derivative
    # The coordinate z runs from the boundary (z=0) to the horizon (z=zh=1)
    def f(z):
        # Handle the edge case at the horizon to prevent sqrt of negative number
        if z >= 1.0:
            return 0.0
        term_under_sqrt = 1.0 - 4.0 * lambda_gb * (1.0 - z**4)
        return (1.0 - np.sqrt(term_under_sqrt)) / (2.0 * lambda_gb)

    def f_prime(z):
        # Derivative of f(z) with respect to z
        if z >= 1.0:
            # The limit as z->1 is well-defined
            return -4.0
        term_under_sqrt = 1.0 - 4.0 * lambda_gb * (1.0 - z**4)
        if term_under_sqrt <= 0:
            return np.nan # Should not be reached in valid range
        return (-4.0 * z**3) / np.sqrt(term_under_sqrt)

    # 3. Define the system of first-order ODEs for the scalar field psi
    # The second-order ODE is:
    # psi'' + (f'/f - 3/z)psi' - (phi^2/f^2 + m^2/(z^2*f))psi = 0
    def ode_system(z, Y, mu):
        psi, d_psi = Y
        
        f_val = f(z)
        # Prevent division by zero very close to the horizon
        if f_val <= 1e-12:
            return [np.nan, np.nan]
            
        fp_val = f_prime(z)
        phi_val = mu * (1.0 - z**2) # Profile for the gauge field A_t

        # Coefficients of the ODE: psi'' + P(z)psi' + Q(z)psi = 0
        P_z = fp_val / f_val - 3.0 / z
        Q_z = - (phi_val**2 / f_val**2) - (m_squared / (z**2 * f_val))

        d2_psi = -P_z * d_psi - Q_z * psi
        return [d_psi, d2_psi]

    # 4. Define a function that shoots from the horizon to the boundary
    # and returns the value of the source term.
    def get_source_term(mu):
        # We want to find mu such that the source term A is zero, where
        # psi(z) ~ A*z + B*z^3 for z -> 0.

        # Integration range starts near the horizon and ends near the boundary
        z_ir = 1.0 - 1e-6
        z_uv = 1e-6

        # Regularity at the horizon sets the initial conditions for psi and its derivative
        psi_ir = 1.0 # The overall scale is arbitrary
        d_psi_ir = -(mu**2) / 4.0 * psi_ir
        Y_ir = [psi_ir, d_psi_ir]
        
        # Integrate the ODE from the IR to the UV
        sol = solve_ivp(
            lambda z, Y: ode_system(z, Y, mu),
            (z_ir, z_uv),
            Y_ir,
            method='RK45',
            dense_output=True,
            atol=1e-9, rtol=1e-9
        )

        # If the solver fails, return a large number to guide the root finder away
        if sol.status != 0:
            return 1e6

        # Extract the solution near the boundary z_uv
        psi_uv, d_psi_uv = sol.sol(z_uv)
        
        # The source term 'A' is found from the asymptotic solution at the boundary
        # A = (3*psi_uv - z_uv*d_psi_uv) / (2*z_uv)
        source = (3 * psi_uv - z_uv * d_psi_uv) / (2 * z_uv)
        return source

    # 5. Use a root-finding algorithm to find mu_c
    # We need to find a search interval [mu_low, mu_high] where get_source_term(mu) changes sign.
    # The value for lambda_gb=0 is ~4.06. For lambda_gb>0, it is expected to be higher.
    try:
        mu_c = brentq(get_source_term, 4.5, 6.0, xtol=1e-6)
        
        # Final result printing
        final_value = mu_c
        equation_string = f"mu_c = {final_value:.4f}"
        
        print(f"For a Gauss-Bonnet coupling of lambda_gb = {lambda_gb}, the operator condenses at a")
        print("critical chemical potential given by the following value:")
        print(equation_string)

    except ValueError:
        print("Could not find the root in the specified interval [4.5, 6.0].")
        print("The calculation might be sensitive or the interval needs adjustment.")
        
# Run the calculation
solve_for_mu_c()
