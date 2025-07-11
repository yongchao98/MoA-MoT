import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation
    in a 5D Einstein-Gauss-Bonnet holographic model.
    """
    # 1. Define Model Parameters
    # Gauss-Bonnet coupling
    lambda_gb = 0.1
    # Mass-squared of the scalar field. For D7-brane instability, this is 0.
    # This corresponds to a dual operator of dimension Delta = 4.
    m_sq = 0.0
    # Charge of the scalar field (can be set to 1 by convention).
    q = 1.0

    # 2. Define Metric Functions
    # The EGB black hole metric function f(z) and its derivative.
    # We set the horizon radius z_h = 1.
    def f(z, lam):
        # A small numerical offset to prevent issues at z=1
        if z > 1.0: z = 1.0
        # The argument of sqrt can become negative due to precision errors near z=1
        arg = 1 - 4 * lam * (1 - z**4)
        if arg < 0: arg = 0
        return (1 / (2 * lam)) * (1 - np.sqrt(arg))

    def f_prime(z, lam):
        if z == 0: return 0
        # A small numerical offset to prevent division by zero at z=1
        arg = 1 - 4 * lam * (1 - z**4)
        if arg <= 0: return 1e-12 # Should return infinity, but sqrt is in denominator
        return (4 * z**3) / np.sqrt(arg)

    # 3. Define the ODE System
    # System for y = [psi, dpsi/dz] from the scalar field's equation of motion.
    def ode_system(z, y, mu, lam, m2):
        psi, dpsi = y
        # Handle z=0 singularity
        if z == 0:
            return [dpsi, 0]

        phi_z = mu * (1 - z**2)
        f_z = f(z, lam)
        if f_z == 0: return [dpsi, 0] # Avoid division by zero at horizon
        f_p_z = f_prime(z, lam)

        # Equation: psi'' + P(z)psi' + Q(z)psi = 0
        p_z = (f_p_z / f_z) - (3 / z)
        q_z = ((q**2 * phi_z**2) / (z**4 * f_z**2)) - (m2 / (z**2 * f_z))
        
        ddpsi = -p_z * dpsi - q_z * psi
        return [dpsi, ddpsi]

    # 4. Define the Shooting Function
    # This function returns the source term of the dual operator, which we want to be zero.
    def get_source_term(mu, lam, m2):
        # Integration range starts near the horizon and ends near the boundary.
        z_start = 1.0 - 1e-6
        z_end = 1e-7

        # Initial conditions for a regular solution at the horizon.
        # psi is normalized to 1, and its derivative is 0.
        y_start = [1.0, 0.0]

        # Solve the ODE using a robust solver.
        sol = solve_ivp(
            ode_system, [z_start, z_end], y_start,
            args=(mu, lam, m2), dense_output=True,
            method='RK45', rtol=1e-10, atol=1e-10
        )

        # Extract solution near the boundary z=0.
        psi_end, dpsi_end = sol.sol(z_end)

        # For m^2=0, solution near boundary is psi(z) = c_0 + c_4*z^4.
        # c_0 is the source, c_4 is the condensate.
        # We find mu that makes the source c_0 zero.
        # c_0 = psi(z) - z*psi'(z)/4 for small z.
        source = psi_end - (z_end * dpsi_end) / 4.0
        return source

    # 5. Find the Root
    # Use a root-finding algorithm (Brent's method) to find mu where the source is zero.
    # We need a bracket [a, b] where the function changes sign.
    # Prior knowledge or simple search shows the root is between 2 and 3.
    try:
        mu_numeric = brentq(get_source_term, 2.0, 3.0, args=(lambda_gb, m_sq))
    except ValueError:
        print("Error: Could not find a root in the given interval.")
        print("The function values at the interval endpoints may not have opposite signs.")
        return

    # 6. Calculate Physical Ratio and Print Results
    # The temperature for z_h=1 is T = 1/pi.
    # The dimensionless ratio is mu_c/T_c = mu_numeric * pi.
    ratio_mu_T = mu_numeric * np.pi
    
    print(f"For a Gauss-Bonnet coupling of lambda_GB = {lambda_gb}:")
    print("The critical chemical potential (mu_c) and critical temperature (T_c) for the condensation are related by the equation:")
    print(f"mu_c = {ratio_mu_T:.4f} * T_c")

# Execute the calculation and print the result.
solve_critical_potential()