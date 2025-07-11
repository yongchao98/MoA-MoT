import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma

def solve_balloon_radius():
    """
    This function solves the given second-order ODE for the radius of a spherical balloon y(t)
    and finds the value of y at t = pi/4.
    """
    
    # 1. Define the initial constant C from y(0)
    # y(0) = (128 * 3^(1/6) * Gamma(2/3))^(-1)
    y0 = (128 * 3**(1/6) * gamma(2/3))**(-1)

    # 2. Define the ODE system for the solver
    # Let Y = [y, y'] = [Y[0], Y[1]]
    # Then Y' = [y', y''] = [Y[1], y'']
    def ode_system(t, Y):
        # Avoid computation at t=0 due to singularity
        if t == 0:
            return [0, 0]

        # Pre-compute trigonometric functions for efficiency and readability
        sin_t = np.sin(t)
        cos_t = np.cos(t)
        tan_t = sin_t / cos_t
        sec_t = 1.0 / cos_t
        t_squared = t * t
        t_cubed = t_squared * t
        t_quad = t_squared * t_squared
        p_t = t_quad + 1 # a common term (t^4+1)

        # Coefficients of the ODE: A2*y'' + A1*y' + A0*y = RHS
        # Original form: 4*(t^4+1)*tan(t)*sec(t)*y'' + ... = (t^4+1)*y*sqrt(sin(t))
        
        # Coefficient of y''(t)
        A2 = 4 * p_t * tan_t * sec_t
        
        # Coefficient of y'(t)
        A1 = 2 * (t_quad + 2 * tan_t * (p_t * tan_t + 8 * t_cubed) + 1) * sec_t
        
        # Coefficient of y(t)
        A0 = 8 * t_squared * (t + 2 * tan_t * (t * tan_t + 3)) * sec_t

        # Right-hand side of the equation
        RHS = p_t * Y[0] * np.sqrt(sin_t)
        
        # Calculate y'' = (RHS - A1*y' - A0*y) / A2
        dY2_dt = (RHS - A1 * Y[1] - A0 * Y[0]) / A2
        
        return [Y[1], dY2_dt]

    # 3. Handle the singularity at t=0
    # Start integration at a small positive time epsilon
    eps = 1e-8
    
    # Use series expansion to find initial values at t = eps
    # y(t) ≈ y0 * (1 + t^(3/2)/6)
    # y'(t) ≈ y0 * (t^(1/2)/4)
    y_eps = y0 * (1 + (1/6) * eps**1.5)
    yp_eps = y0 * (1/4) * eps**0.5
    
    initial_conditions_eps = [y_eps, yp_eps]
    
    # 4. Solve the ODE numerically
    t_span = [eps, np.pi/4]
    # 'Radau' is a good method for stiff problems which can arise near singularities
    sol = solve_ivp(ode_system, t_span, initial_conditions_eps, method='Radau', dense_output=True)
    
    # 5. Extract and print the final radius at t = pi/4
    radius_at_pi_over_4 = sol.y[0, -1]
    
    print(f"The initial radius y(0) is: {y0}")
    print(f"The radius of the balloon at t = pi/4 is y(pi/4) = {radius_at_pi_over_4}")

solve_balloon_radius()