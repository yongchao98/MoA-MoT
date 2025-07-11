import numpy as np
from scipy.integrate import solve_ivp
from math import gamma, pi

def ode_system(t, z):
    """
    Defines the system of first-order ODEs.
    z = [y, y']
    """
    y, y_p = z
    
    # Handle the edge case t=0 if the solver ever tries it.
    if t == 0:
        return [0, 0]

    # Pre-calculate trigonometric terms for efficiency and clarity
    cos_t = np.cos(t)
    sin_t = np.sin(t)
    tan_t = sin_t / cos_t
    sec_t = 1 / cos_t
    
    t_sq = t**2
    t_cub = t**3
    t_quad = t**4
    t4_plus_1 = t_quad + 1

    # Define the coefficients of the ODE: A*y'' + B*y' + C*y = D*y
    A = 4 * t4_plus_1 * tan_t * sec_t
    
    B = 2 * sec_t * (t_quad + 2 * tan_t * (t4_plus_1 * tan_t + 8 * t_cub) + 1)
    
    C = 8 * t_sq * sec_t * (t + 2 * tan_t * (t * tan_t + 3))
    
    # The term on the right-hand side
    D = t4_plus_1 * np.sqrt(sin_t)
    
    # Calculate y'' = (1/A) * ((D-C)*y - B*y')
    y_pp = (1/A) * ((D - C) * y - B * y_p)
    
    return [y_p, y_pp]

# Given initial value y(0)
y0_val = (128 * 3**(1/6) * gamma(2/3))**(-1)

# We start the integration at a very small time t0 > 0 to avoid the singularity at t=0.
t0 = 1e-8

# From series expansion around t=0, y(t) approx y(0)*(1 + t^(3/2)/9)
# So we can approximate y(t0) and y'(t0).
y_t0 = y0_val * (1 + (1/9) * t0**1.5)
yp_t0 = y0_val * (1/6) * t0**0.5

# Set the initial conditions for the solver
z0 = [y_t0, yp_t0]

# Define the time span for the integration
t_span = [t0, pi/4]

# Solve the ODE
# We use a high precision solver and tight tolerances for an accurate result.
sol = solve_ivp(ode_system, t_span, z0, method='RK45', rtol=1e-12, atol=1e-12)

# The result y(pi/4) is the last value of the first component of the solution vector.
y_pi_over_4 = sol.y[0, -1]

# The final equation is y(pi/4) = value. We print the value.
print(y_pi_over_4)