import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma

# Define the coefficients of the ODE Ay'' + By' + Cy = Dy
def A(t):
    if t == 0:
        return 0
    return 4 * (t**4 + 1) * np.tan(t) * (1/np.cos(t))

def B(t):
    tan_t = np.tan(t)
    sec_t = 1/np.cos(t)
    term1 = t**4 + 1
    term2 = 2 * tan_t * ((t**4 + 1) * tan_t + 8 * t**3)
    return 2 * (term1 + term2) * sec_t

def C(t):
    tan_t = np.tan(t)
    sec_t = 1/np.cos(t)
    return 8 * t**2 * (t + 2 * tan_t * (t * tan_t + 3)) * sec_t

def D(t):
    return (t**4 + 1) * np.sqrt(np.sin(t))

# Define the ODE system for the solver, Y' = F(t, Y), where Y = [y, y']
def ode_system(t, Y):
    y, y_prime = Y
    if t == 0:
        # At t=0, A(t) is 0. This is a singular point.
        # From asymptotic analysis 4ty'' + 2y' ~ y sqrt(t)
        # As t->0+, y'' ~ (y sqrt(t) - 2y') / 4t
        # Using y'(t) ~ y(0)/4 * t^(1/2), y''(t) ~ y(0)/8 * t^(-1/2) -> diverges
        # The solver will not start at t=0, so this case is theoretical.
        # In practice we start from a small epsilon.
        return [0, 0] # Placeholder, should not be called

    y_double_prime = (D(t) * y - B(t) * y_prime - C(t) * y) / A(t)
    return [y_prime, y_double_prime]

# Initial conditions
y0_const = 128 * 3**(1/6) * gamma(2/3)
y0 = 1 / y0_const
yp0 = 0

# The ODE is singular at t=0. Start integration from a small positive t_start.
t_start = 1e-8

# Use asymptotic solution y(t) ~= y0 * (1 + 1/6 * t^(3/2)) for small t to find initial conditions at t_start.
y_start = y0 * (1 + (1/6) * t_start**1.5)
yp_start = y0 * (1/4) * t_start**0.5

# Time span for integration
t_end = np.pi / 4
t_span = [t_start, t_end]

# Solve the ODE
sol = solve_ivp(ode_system, t_span, [y_start, yp_start], dense_output=True, rtol=1e-8, atol=1e-10)

# Get the solution at t=pi/4
t_eval = np.pi / 4
y_eval, yp_eval = sol.sol(t_eval)

# Calculate y'' at t=pi/4 from the ODE itself to verify
_, ypp_eval = ode_system(t_eval, [y_eval, yp_eval])

# Get coefficients at t=pi/4
A_eval = A(t_eval)
B_eval = B(t_eval)
C_eval = C(t_eval)
D_eval = D(t_eval)

# Display the equation with the computed values
print("The differential equation at t = pi/4:")
lhs = A_eval * ypp_eval + B_eval * yp_eval + C_eval * y_eval
rhs = D_eval * y_eval
print(f"{A_eval:.5f} * {ypp_eval:.5f} + {B_eval:.5f} * {yp_eval:.5f} + {C_eval:.5f} * {y_eval:.5f} = {D_eval:.5f} * {y_eval:.5f}")
print(f"{lhs:.5f} = {rhs:.5f}")

# Print the final result
print(f"\nThe radius y(t) at t=pi/4 is: {y_eval}")

print(f"<<<{y_eval}>>>")