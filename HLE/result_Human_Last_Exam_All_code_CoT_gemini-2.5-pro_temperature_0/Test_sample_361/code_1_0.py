import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.misc import derivative

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a 5D
    Einstein-Gauss-Bonnet holographic model.
    """
    # --- 1. Define Parameters ---
    LGB = 0.1  # Gauss-Bonnet coupling
    RH = 1.0   # Horizon radius (can be set to 1 by scaling)
    M2 = -3.0  # Mass squared of the scalar field for D3/D7
    Q = 1.0    # Charge of the scalar field

    # Integration parameters
    R_START = RH + 1e-5
    R_MAX = 200.0

    # --- 2. Define Background Functions ---
    def f(r, lgb=LGB):
        """EGB black hole metric function f(r)."""
        if r < RH:
            return 0
        # Avoid sqrt of negative for r slightly less than RH due to precision
        arg_sqrt = 1 - 4 * lgb * (1 - (RH**4 / r**4))
        if arg_sqrt < 0:
            arg_sqrt = 0
        return (r**2 / (2 * lgb)) * (1 - np.sqrt(arg_sqrt))

    def f_prime(r, lgb=LGB):
        """Numerical derivative of f(r)."""
        return derivative(lambda x: f(x, lgb), r, dx=1e-6)

    # --- 3. Define the ODE System ---
    def ode_system(r, y, mu):
        """
        Defines the system of ODEs for y = [psi, psi_prime].
        The EOM is: psi'' + (3/r + f'/f)psi' - (q^2*phi^2/f^2 + m^2/f)psi = 0
        """
        psi, psi_p = y
        fr = f(r)
        if fr <= 1e-12:  # At or very near the horizon
            return [0.0, 0.0]
        
        fr_p = f_prime(r)
        phi = mu * (1 - RH**2 / r**2)
        
        # Rearrange EOM to solve for psi''
        psi_pp = - (3/r + fr_p/fr) * psi_p + (Q**2 * phi**2 / fr**2 + M2 / fr) * psi
        
        return [psi_p, psi_pp]

    # --- 4. Setup the Shooting Method ---
    # Calculate asymptotic exponents to check the boundary condition
    # As r -> infinity, f(r) -> c2 * r^2
    c2 = (1 / (2 * LGB)) * (1 - np.sqrt(1 - 4 * LGB))
    m2_eff = M2 / c2
    # The field behaves as psi ~ A*r^delta_plus + B*r^delta_minus
    # delta are roots of Delta^2 + 4*Delta - m2_eff = 0
    delta_plus = -2 + np.sqrt(4 + m2_eff)   # Slower fall-off (source)
    delta_minus = -2 - np.sqrt(4 + m2_eff)  # Faster fall-off (VEV)

    def shoot(mu):
        """
        Objective function for the root finder. It returns the coefficient
        of the source term, which should be zero for condensation.
        """
        # Initial conditions corresponding to a regular solution at the horizon
        y0 = [1.0, 0.0]
        
        # Solve the ODE from near the horizon to the boundary
        sol = solve_ivp(
            lambda r, y: ode_system(r, y, mu),
            [R_START, R_MAX],
            y0,
            method='RK45',
            dense_output=True,
            atol=1e-9, rtol=1e-9
        )
        
        # Extract solution at the boundary (R_MAX)
        r_f = R_MAX
        psi_f, psip_f = sol.sol(r_f)
        
        # We want the coefficient of the source term (r^delta_plus) to be zero.
        # The numerator of this coefficient is proportional to:
        source_coeff_numerator = r_f * psip_f - delta_minus * psi_f
        return source_coeff_numerator

    # --- 5. Find the Root ---
    try:
        # Search for the lowest mu > 0 where the source coefficient is zero
        mu_c = brentq(shoot, 1.0, 8.0, xtol=1e-6)
        
        # Print the final result in a full sentence.
        # The problem asks to output each number in the final equation.
        # We interpret this as clearly stating all parameters and the result.
        print(f"For a D3/D7 holographic model with Gauss-Bonnet coupling lambda_GB = {LGB} and scalar field mass squared m^2 = {M2}, the critical chemical potential is:")
        print(f"mu_c = {mu_c}")

    except ValueError:
        print("Root not found in the given interval. The shooting method failed.")

if __name__ == '__main__':
    solve_critical_potential()