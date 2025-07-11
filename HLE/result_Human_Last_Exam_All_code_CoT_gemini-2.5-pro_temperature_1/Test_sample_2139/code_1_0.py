import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma

def solve_balloon_radius():
    """
    Solves the given differential equation for the radius of a spherical balloon y(t)
    and finds its value at t = pi/4.
    """

    # Define the system of first-order ODEs from the given second-order ODE
    # y' = z
    # z' = y'' = ( (t^4+1)sqrt(sin(t))y - P1(t)z - P2(t)y ) / P0(t)
    def ode_system(t, Y):
        y, z = Y # Y = [y, y']

        # The coefficients are singular or zero at t=0, so this function is only valid for t > 0
        if t <= 0:
            return [0, 0]

        tan_t = np.tan(t)
        sec_t = 1.0 / np.cos(t)
        
        # Coefficient of y''(t)
        P0 = 4 * (t**4 + 1) * tan_t * sec_t

        # Coefficient of y'(t)
        P1 = 2 * (t**4 + 2 * tan_t * ((t**4 + 1) * tan_t + 8 * t**3) + 1) * sec_t

        # Coefficient of y(t) on the LHS
        P2 = 8 * t**2 * (t + 2 * tan_t * (t * tan_t + 3)) * sec_t

        # Term with y(t) on the RHS
        RHS_term = (t**4 + 1) * np.sqrt(np.sin(t)) * y
        
        # Calculate y'' (which is z')
        y_double_prime = (RHS_term - P1 * z - P2 * y) / P0
        
        return [z, y_double_prime]

    # Calculate the initial condition y(0)
    y0_val = (128 * 3**(1/6) * gamma(2/3))**(-1)

    # The ODE is singular at t=0. We start the integration from a small positive value t_start.
    t_start = 1e-9

    # Based on asymptotic analysis near t=0, we find the initial values for the solver.
    # y(t_start) is approximately y(0)
    # y'(t_start) is approximately y(0)/3 * sqrt(t_start)
    y_start = y0_val
    y_prime_start = y0_val / 3.0 * np.sqrt(t_start)
    
    # Initial state vector at t_start
    Y_start = [y_start, y_prime_start]

    # Time span for the integration
    t_end = np.pi / 4
    t_span = [t_start, t_end]

    # Solve the ODE using a robust method (RK45) with tight tolerances
    sol = solve_ivp(
        ode_system, 
        t_span, 
        Y_start, 
        method='RK45', 
        dense_output=True,
        atol=1e-10, 
        rtol=1e-10
    )

    # Extract the solution at t = pi/4
    y_at_pi_over_4 = sol.sol(t_end)[0]
    
    # Print the result
    # The format "output each number in the final equation" is interpreted as
    # printing the final calculated value in a descriptive equation format.
    print(f"y({t_end:.4f}) = {y_at_pi_over_4}")

solve_balloon_radius()