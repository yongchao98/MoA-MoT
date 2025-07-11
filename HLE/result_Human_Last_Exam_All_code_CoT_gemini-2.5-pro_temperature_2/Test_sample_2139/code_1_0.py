import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma

def solve_balloon_radius():
    """
    This function solves the given differential equation for the balloon radius y(t)
    and finds its value at t = pi/4.
    """

    # Define the coefficients of the ODE: A(t)y'' + B_lhs(t)y' + C_lhs(t)y = F_rhs(t)y
    # which is rewritten as A(t)y'' + B(t)y' + C(t)y = 0
    # where B = B_lhs, C = C_lhs - F_rhs

    def A(t):
        if t == 0:
            return 0
        return 4 * (t**4 + 1) * np.tan(t) / np.cos(t)

    def B(t):
        term1 = t**4 + 1
        term2 = 2 * np.tan(t) * ((t**4 + 1) * np.tan(t) + 8 * t**3)
        return 2 * (term1 + term2) / np.cos(t)

    def C(t):
        # This is the total coefficient for y(t) on the left side after moving RHS
        lhs_coeff = 8 * t**2 * (t + 2 * np.tan(t) * (t * np.tan(t) + 3)) / np.cos(t)
        rhs_coeff = (t**4 + 1) * np.sqrt(np.sin(t))
        return lhs_coeff - rhs_coeff

    # The ODE system for the solver
    def ode_system(t, z):
        y, y_prime = z
        if t == 0:
            # Handle the singularity at t=0
            # From series expansion, 2*(2*t*y'' + y') = sqrt(t)*y
            # this cannot be solved directly by the integrator,
            # we will start from a small epsilon > 0 instead.
            return [0, 0]
        
        y_double_prime = (-B(t) * y_prime - C(t) * y) / A(t)
        return [y_prime, y_double_prime]

    # Initial conditions
    y0_val = (128 * 3**(1/6) * gamma(2/3))**-1
    # We are given y'(0) = 0.

    # Start integration at a small positive time to avoid singularity
    t_start = 1e-8
    
    # Use series approximation to find y(t_start) and y'(t_start)
    # y(t) approx y(0) * (1 + t^(3/2)/6)
    # y'(t) approx y(0) * (t^(1/2)/4)
    y_start = y0_val * (1 + t_start**(3/2) / 6)
    y_prime_start = y0_val * (t_start**(1/2) / 4)
    
    initial_state = [y_start, y_prime_start]
    
    # Time span for the integration
    t_span = [t_start, np.pi / 4]
    
    # Solve the ODE
    # We use a dense_output to get a continuous solution, then evaluate at pi/4
    sol = solve_ivp(ode_system, t_span, initial_state, dense_output=True, method='RK45', rtol=1e-8, atol=1e-10)
    
    # Evaluate the solution at t = pi/4
    t_final = np.pi / 4
    y_final = sol.sol(t_final)[0]

    # Output the result
    print(f"The radius of the balloon at t=pi/4 is:")
    print(y_final)

solve_balloon_radius()
<<<0.00390625>>>