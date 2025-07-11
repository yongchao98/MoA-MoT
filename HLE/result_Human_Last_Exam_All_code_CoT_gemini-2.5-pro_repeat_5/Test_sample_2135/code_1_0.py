import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the differential equation dy/dt = f(t, y)
# The function takes t (time) and y (angle) and returns the rate of change dy/dt.
def pendulum_ode(t, y):
    """
    Represents the differential equation dy/dt = sin(y(t) - t) - 4.
    y is expected to be an array-like object from solve_ivp, so we use y[0].
    """
    return np.sin(y[0] - t) - 4

# Step 2: Set up the problem parameters
# Initial condition: y(0) = 0
y0 = [0]
# Target time: t = pi/6
t_target = np.pi / 6
# Time interval for the solver
t_span = [0, t_target]

# Step 3: Solve the differential equation numerically to find y(pi/6)
# We use solve_ivp from the SciPy library for an accurate solution.
# We ask the solver to return the solution only at our target time t_target.
sol = solve_ivp(pendulum_ode, t_span, y0, t_eval=[t_target])

# Extract the calculated value of y at t = pi/6
# The result is in sol.y, which is a 2D array. For our single equation and single
# evaluation point, the value is at sol.y[0, 0].
y_at_target = sol.y[0, 0]

# Step 4: Calculate the instantaneous rate of change (dy/dt) at t = pi/6
# We plug the values of t and y(t) back into the original differential equation.
# The rate of change is dy/dt = sin(y(pi/6) - pi/6) - 4
rate_of_change = np.sin(y_at_target - t_target) - 4

# Step 5: Print the results in a clear, step-by-step format
print("The problem is to find the instantaneous rate of change dy/dt at t = π/6.")
print("The motion is described by the differential equation: dy/dt = sin(y(t) - t) - 4")
print(f"The initial condition is y(0) = {y0[0]}.")
print("\nFirst, we need to find the value of y(t) at t = π/6.")
print("We solve the differential equation numerically from t=0 to t=π/6.")
print(f"\nAt t = π/6 ≈ {t_target:.6f}, the calculated value of the angle is y(π/6) ≈ {y_at_target:.6f}.")
print("\nNow, we can calculate the instantaneous rate of change dy/dt at this point.")
print("The formula is: dy/dt = sin(y(π/6) - π/6) - 4")
print(f"Substituting the values: dy/dt = sin({y_at_target:.6f} - {t_target:.6f}) - 4")
argument_of_sin = y_at_target - t_target
print(f"dy/dt = sin({argument_of_sin:.6f}) - 4")
value_of_sin = np.sin(argument_of_sin)
print(f"dy/dt = {value_of_sin:.6f} - 4")
print(f"\nThe final result is:")
print(f"dy/dt = {rate_of_change:.6f}")