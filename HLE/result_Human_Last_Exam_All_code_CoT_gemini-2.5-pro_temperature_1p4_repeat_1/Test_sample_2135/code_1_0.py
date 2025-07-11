import numpy as np
from scipy.integrate import solve_ivp

# Define the differential equation dy/dt = f(t, y)
def pendulum_motion(t, y):
    """
    Defines the differential equation for the pendulum's motion.
    The input y is a list or array, so we access its first element y[0].
    """
    return np.sin(y[0] - t) - 4

# Set up the problem parameters
# Initial condition: y(0) = 0
y0 = [0]
# The time at which we need to calculate the rate of change
t_target = np.pi / 6
# The time span for the ODE solver [start_time, end_time]
t_span = [0, t_target]

# Step 1: Numerically solve the ODE to find the value of y at t = pi/6.
# We instruct the solver to return the solution specifically at our target time.
solution = solve_ivp(pendulum_motion, t_span, y0, t_eval=[t_target])

# Extract the value of y(pi/6) from the solution object.
y_at_target_t = solution.y[0, 0]

# Step 2: Calculate the instantaneous rate of change dy/dt at t = pi/6
# using the value of y(pi/6) we just found.
rate_of_change = pendulum_motion(t_target, [y_at_target_t])

# Print the step-by-step calculation
print("The differential equation is: dy/dt = sin(y(t) - t) - 4")
print(f"The initial condition is y(0) = {y0[0]}.")
print(f"We need to find dy/dt at t = π/6 ≈ {t_target:.6f}.")
print("\nFirst, we solve for y(π/6) by integrating the differential equation:")
print(f"The numerical solution gives y(π/6) ≈ {y_at_target_t:.6f}.")

print("\nNow, we calculate the instantaneous rate of change using the formula:")
print("dy/dt = sin(y(π/6) - π/6) - 4")
print("\nPlugging in the calculated values:")
# For clarity, let's calculate the terms separately
argument_of_sin = y_at_target_t - t_target
sin_term = np.sin(argument_of_sin)

print(f"dy/dt = sin({y_at_target_t:.6f} - {t_target:.6f}) - 4")
print(f"dy/dt = sin({argument_of_sin:.6f}) - 4")
print(f"dy/dt = {sin_term:.6f} - 4")
print(f"dy/dt = {rate_of_change:.6f}")