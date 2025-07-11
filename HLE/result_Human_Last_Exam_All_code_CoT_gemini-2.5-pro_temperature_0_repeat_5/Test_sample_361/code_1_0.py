import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using the shooting method.
    """
    # 1. Define constants and parameters
    LAM_GB = 0.1  # Gauss-Bonnet coupling
    R_H = 1.0     # Horizon radius (set to 1 to fix the energy scale)

    # 2. Define the metric and gauge field functions
    def f_metric(r, lam, rh):
        """
        The EGB black hole metric function f(r).
        """
        # The term inside the square root must be non-negative.
        # For r >= rh, this is guaranteed if 1 - 4*lam >= 0.
        # lam = 0.1 is safe.
        term_in_sqrt = 1 - 4 * lam * (1 - (rh**4 / r**4))
        return (r**2 / (2 * lam)) * (1 - np.sqrt(term_in_sqrt))

    def f_prime_metric(r, lam, rh):
        """
        Numerical derivative of the metric function f(r).
        """
        h = 1e-6
        return (f_metric(r + h, lam, rh) - f_metric(r - h, lam, rh)) / (2 * h)

    # 3. Define the ODE system from the Equation of Motion
    def ode_system(r, y, mu, lam, rh):
        """
        Represents the second-order ODE for the scalar field as a system of
        first-order ODEs. y = [phi, dphi/dr].
        """
        phi, dphi = y
        
        f_val = f_metric(r, lam, rh)
        fp_val = f_prime_metric(r, lam, rh)
        
        # To avoid division by zero right at the horizon
        if abs(f_val) < 1e-12:
            f_val = 1e-12

        # Coefficients of the ODE: phi'' + P(r)phi' + Q(r)phi = 0
        P_r = fp_val / f_val + 3 / r
        
        # A_t(r) = mu * (1 - rh^2/r^2)
        At_val = mu * (1 - rh**2 / r**2)
        Q_r = At_val**2 / f_val**2
        
        d2phi = -P_r * dphi - Q_r * phi
        return [dphi, d2phi]

    # 4. Define the objective function for the shooting method
    def objective_function(mu):
        """
        This function performs the "shot" and returns a value that is zero
        when the boundary condition at infinity is met.
        """
        # Integration parameters
        r_start = R_H + 1e-5
        r_end = 800.0  # A sufficiently large radius for the boundary

        # Initial conditions at the horizon, derived from regularity
        phi_0 = 1.0  # The overall scale is arbitrary
        # phi'(rh) = - (mu^2 * (1 - 4*lam) / (4*rh^4)) * phi(rh)
        dphi_0 = - (mu**2 * (1 - 4 * LAM_GB) / (4 * R_H**4)) * phi_0
        y0 = [phi_0, dphi_0]

        # Numerically integrate the ODE
        sol = solve_ivp(
            ode_system,
            [r_start, r_end],
            y0,
            args=(mu, LAM_GB, R_H),
            method='RK45',
            dense_output=True
        )

        # Extract the solution at the boundary r_end
        phi_end, dphi_end = sol.sol(r_end)

        # The asymptotic solution is phi(r) ~ c_source + c_condensate/r^4.
        # We want the solution with zero source, c_source = 0.
        # c_source can be extracted from the solution via:
        # c_source = phi(r) + r * phi'(r) / 4
        c_source = phi_end + r_end * dphi_end / 4.0
        
        return c_source

    # 5. Find the root of the objective function
    # Based on literature, the value for lam=0.1 is expected to be around 8.3.
    # We choose a bracket around this value for the root finder.
    search_bracket = [7.0, 10.0]
    try:
        result = root_scalar(objective_function, bracket=search_bracket, method='brentq')
        mu_c = result.root
        print("The equation for the critical chemical potential is:")
        print(f"μ_c = {mu_c:.4f}")
    except ValueError:
        # This error occurs if the objective function has the same sign at both
        # ends of the bracket, meaning no root was found.
        val_a = objective_function(search_bracket[0])
        val_b = objective_function(search_bracket[1])
        print("Could not find the root in the specified bracket.")
        print(f"Objective function at μ={search_bracket[0]}: {val_a:.4f}")
        print(f"Objective function at μ={search_bracket[1]}: {val_b:.4f}")
        print("Please try adjusting the search_bracket.")

if __name__ == '__main__':
    solve_critical_potential()