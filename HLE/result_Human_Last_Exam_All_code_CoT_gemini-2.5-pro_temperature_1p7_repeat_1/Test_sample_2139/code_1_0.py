import numpy as np
from scipy.special import gamma
from scipy.integrate import solve_ivp

def ode_system(t, z):
    """
    Defines the system of first-order ODEs.
    z[0] = y
    z[1] = y'
    Returns z' = [y', y'']
    """
    y, yp = z[0], z[1]

    # Handle t=0 case to avoid division by zero, although the solver won't call it at t=0.
    if t == 0:
        return [0, 0] # This should not be reached with the chosen integration interval

    # Coefficients of the ODE: A*y'' + B*y' + C*y = G*y
    # A(t)
    coeff_ypp = 4 * (t**4 + 1) * np.tan(t) * (1/np.cos(t))
    
    # B(t)
    coeff_yp = 2 * (t**4 + 2 * np.tan(t) * ((t**4 + 1) * np.tan(t) + 8 * t**3) + 1) * (1/np.cos(t))

    # C(t) (from y term on LHS)
    coeff_y_lhs = 8 * t**2 * (t + 2 * np.tan(t) * (t * np.tan(t) + 3)) * (1/np.cos(t))

    # G(t) (from y term on RHS)
    coeff_y_rhs = (t**4 + 1) * np.sqrt(np.sin(t))
    
    # Calculate y'' = (G*y - B*y' - C*y) / A
    ypp = (coeff_y_rhs * y - coeff_yp * yp - coeff_y_lhs * y) / coeff_ypp

    return [yp, ypp]

def solve_balloon_radius():
    """
    Solves for the radius of the spherical balloon y(t) at t=pi/4.
    """
    # Calculate the initial value y(0)
    y0 = 1 / (128 * (3**(1/6)) * gamma(2/3))

    # We cannot start at t=0 because the equation is singular.
    # We use an asymptotic expansion for small t to find initial conditions at t_start.
    t_start = 1e-8

    # y(t) approx y(0) * (1 + t^(3/2)/6) for small t
    y_start = y0 * (1 + (t_start**1.5) / 6)
    
    # y'(t) approx y(0) * t^(1/2) / 4 for small t
    yp_start = y0 * (t_start**0.5) / 4
    
    z0 = [y_start, yp_start]

    # Define the time span for the integration
    t_end = np.pi / 4
    t_span = [t_start, t_end]
    
    # Use solve_ivp to find the solution
    sol = solve_ivp(ode_system, t_span, z0, dense_output=True, rtol=1e-8, atol=1e-8)

    # Get the value at t=pi/4
    y_pi_over_4 = sol.sol(t_end)[0]
    
    # Output the result as an equation
    print(f"y(pi/4) = {y_pi_over_4}")

# Execute the solver
solve_balloon_radius()