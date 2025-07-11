import numpy as np
from scipy.integrate import solve_ivp

def pendulum_motion(t, y):
  """
  Defines the differential equation for the pendulum's motion:
  dy/dt = sin(y(t) - t) - 4
  """
  return np.sin(y - t) - 4

# The time 't' at which we want to find the rate of change
t_eval_point = np.pi / 6

# The initial condition: y(0) = 0
y0 = [0]

# The time interval for the solver: from t=0 to t=pi/6
t_span = [0, t_eval_point]

# Step 1: Solve the initial value problem to find y(π/6)
# We instruct the solver to find the solution specifically at t_eval_point.
solution = solve_ivp(
    fun=pendulum_motion,
    t_span=t_span,
    y0=y0,
    t_eval=[t_eval_point]
)

# Extract the calculated value of y at t = π/6
y_at_t_eval = solution.y[0, 0]

# Step 2: Calculate the rate of change dy/dt at t = π/6
# by substituting the found value of y(π/6) into the equation.
rate_of_change = pendulum_motion(t_eval_point, y_at_t_eval)

# --- Output the results step-by-step ---
print(f"The problem is to find dy/dt at t = π/6 for the equation: dy/dt = sin(y(t) - t) - 4, with y(0) = 0.")
print("\nFirst, we solve for y(t) at t = π/6.")
print(f"Using a numerical solver, we find that at t = {t_eval_point:.6f}, the value of y(t) is {y_at_t_eval:.6f}.")

print("\nNow we substitute these values into the equation for the rate of change:")
# Print the equation with the numerical values
print(f"dy/dt = sin(y(π/6) - π/6) - 4")
print(f"dy/dt = sin({y_at_t_eval:.6f} - {t_eval_point:.6f}) - 4")

# Calculate and print the intermediate steps
argument_of_sin = y_at_t_eval - t_eval_point
value_of_sin = np.sin(argument_of_sin)

print(f"dy/dt = sin({argument_of_sin:.6f}) - 4")
print(f"dy/dt = {value_of_sin:.6f} - 4")

# Print the final result
print("\nThe final calculated instantaneous rate of change is:")
print(f"dy/dt = {rate_of_change:.6f}")