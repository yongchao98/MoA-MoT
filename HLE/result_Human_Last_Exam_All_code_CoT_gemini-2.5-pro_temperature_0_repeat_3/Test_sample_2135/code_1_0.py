import numpy as np
from scipy.integrate import solve_ivp

def pendulum_motion(t, y):
    """
    Defines the differential equation for the pendulum's motion.
    dy/dt = sin(y(t) - t) - 4
    """
    return np.sin(y - t) - 4

# Set the time at which we want to find the rate of change
t_eval_point = np.pi / 6

# Set the initial condition: y(0) = 0
initial_condition = [0]

# Define the time span for the integration from t=0 to t=pi/6
time_span = [0, t_eval_point]

# Numerically solve the ODE to find the value of y at t = pi/6
# The 't_eval' argument ensures we get the solution at our specific point.
solution = solve_ivp(
    fun=pendulum_motion,
    t_span=time_span,
    y0=initial_condition,
    t_eval=[t_eval_point]
)

# Extract the calculated value of y(pi/6) from the solution object
y_at_pi_over_6 = solution.y[0, 0]

# Now, calculate the instantaneous rate of change (dy/dt) at t = pi/6
# by plugging the values back into the original equation.
rate_of_change = pendulum_motion(t_eval_point, y_at_pi_over_6)

# Print the results, showing each number in the final equation
print("To find the instantaneous rate of change dy/dt at t = pi/6, we use the equation:")
print("dy/dt = sin(y(t) - t) - 4\n")

print(f"First, we solve for y(t) at t = pi/6 (approx. {t_eval_point:.5f}).")
print(f"The numerical solution gives y(pi/6) = {y_at_pi_over_6:.5f}\n")

print("Next, we substitute these values into the equation for dy/dt:")
print(f"dy/dt = sin({y_at_pi_over_6:.5f} - {t_eval_point:.5f}) - 4")
print(f"dy/dt = sin({y_at_pi_over_6 - t_eval_point:.5f}) - 4")
print(f"dy/dt = {np.sin(y_at_pi_over_6 - t_eval_point):.5f} - 4")
print(f"dy/dt = {rate_of_change:.5f}\n")

print(f"The instantaneous rate of change at t = pi/6 is {rate_of_change:.5f}.")