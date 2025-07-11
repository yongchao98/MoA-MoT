import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma

def solve_balloon_radius():
    """
    Solves the given ODE for the radius of a spherical balloon and finds its value at t=pi/4.
    The solution is found by numerically integrating the differential equation using scipy.integrate.solve_ivp.
    Due to a singularity at t=0, the integration starts from a small epsilon, with initial conditions
    derived from a small-t expansion of the ODE.
    """

    # Calculate the initial condition y(0) from the given expression.
    # y(0) = C0 = (128 * 3^(1/6) * Gamma(2/3))^(-1)
    C0 = 1 / (128 * 3**(1/6) * gamma(2/3))

    # Define the system of first-order ODEs for solve_ivp.
    # z = [y, y'] -> z' = [y', y'']
    def ode_system(t, z):
        y, yp = z
        
        # This check is to avoid division by zero if t is exactly 0.
        if t == 0:
            # The ODE is singular at t=0, so the solver should not evaluate here.
            # Returning nan ensures any accidental call is flagged.
            return [np.nan, np.nan]

        sin_t = np.sin(t)
        cos_t = np.cos(t)
        tan_t = sin_t / cos_t
        sec_t = 1 / cos_t
        t2 = t * t
        t3 = t2 * t
        t4 = t3 * t

        # Rewrite the ODE as y'' = (F(t, y) - Q(t)y' - R(t)y) / P(t)
        # where the original equation is P(t)y'' + Q_orig(t)y' + R_orig(t)y = F_orig(t)y
        # y''(t) = [ (t^4+1)y(t)sqrt(sin(t)) - Q_orig(t)y' - R_orig(t)y ] / P(t)

        # Coefficient of y''(t)
        P_t = 4 * (t4 + 1) * tan_t * sec_t
        
        # Coefficient of y'(t)
        Q_orig_t = 2 * (t4 + 2 * tan_t * ((t4 + 1) * tan_t + 8 * t3) + 1) * sec_t
        
        # Coefficient of y(t) on the LHS
        R_orig_t = 8 * t2 * (t + 2 * tan_t * (t * tan_t + 3)) * sec_t
        
        # Term on the RHS
        F_t = (t4 + 1) * y * np.sqrt(sin_t)
        
        # Calculate y''(t)
        ypp = (F_t - Q_orig_t * yp - R_orig_t * y) / P_t
        
        return [yp, ypp]

    # To handle the singularity at t=0, start integration from a small epsilon.
    # We find y(epsilon) and y'(epsilon) using a small-t approximation.
    # For small t: y(t) ≈ C0 * (1 + t^(3/2)/6) and y'(t) ≈ C0 * t^(1/2)/4.
    epsilon = 1e-9
    
    y_eps = C0 * (1 + epsilon**(1.5) / 6.0)
    yp_eps = C0 * (epsilon**(0.5) / 4.0)
    
    z0 = [y_eps, yp_eps]

    # Define the time interval for the solver.
    t_start = epsilon
    t_end = np.pi / 4
    t_span = [t_start, t_end]

    # We want the solution only at the endpoint t=pi/4.
    t_eval = [t_end]
    
    # Solve the Initial Value Problem.
    sol = solve_ivp(
        ode_system, 
        t_span, 
        z0, 
        method='LSODA', 
        t_eval=t_eval,
        rtol=1e-9,  # Relative tolerance
        atol=1e-11  # Absolute tolerance
    )

    # Extract the final result for y(pi/4).
    y_pi_over_4 = sol.y[0, -1]

    # Print the final equation and its numerical components as requested.
    # The final equation is "y(pi/4) = value".
    print(f"y({np.pi/4}) = {y_pi_over_4}")

solve_balloon_radius()