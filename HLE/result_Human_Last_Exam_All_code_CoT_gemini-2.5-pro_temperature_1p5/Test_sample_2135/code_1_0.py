import numpy as np
from scipy.integrate import solve_ivp

# Define the differential equation dy/dt = f(t, y)
def pendulum_ode(t, y):
  """
  Defines the differential equation for the pendulum's motion.
  y is expected to be a list or array, so we access its first element.
  """
  return np.sin(y[0] - t) - 4

# Set up the parameters for the problem
t_eval = np.pi / 6  # The time at which to evaluate the rate of change
t_span = [0, t_eval] # The time interval for the solver
y0 = [0.0]           # The initial condition y(0) = 0

# 1. Numerically solve the ODE to find y(π/6)
# We ask the solver to find the solution specifically at t_eval.
solution = solve_ivp(pendulum_ode, t_span, y0, t_eval=[t_eval], dense_output=True)
y_at_pi_6 = solution.y[0, 0]

# 2. Calculate the instantaneous rate of change at t = π/6
# dy/dt = sin(y(π/6) - π/6) - 4
rate_of_change = pendulum_ode(t_eval, [y_at_pi_6])

# 3. Display the results and the final equation
print("Step 1: Find the value of y(t) at t = π/6 by solving the ODE.")
print(f"The time t is: {t_eval:.6f} radians (π/6)")
print(f"The calculated angle y(t) at this time is: {y_at_pi_6:.6f} radians")
print("\nStep 2: Calculate the instantaneous rate of change dy/dt at t = π/6.")
print("The formula is: dy/dt = sin(y(t) - t) - 4")
print("\nSubstituting the values into the equation:")

arg = y_at_pi_6 - t_eval
sin_arg = np.sin(arg)

print(f"dy/dt = sin({y_at_pi_6:.6f} - {t_eval:.6f}) - 4")
print(f"dy/dt = sin({arg:.6f}) - 4")
print(f"dy/dt = {sin_arg:.6f} - 4")
print(f"dy/dt = {rate_of_change:.6f}")

print("\nFinal equation with all numbers:")
print(f"sin({y_at_pi_6:.6f} - {t_eval:.6f}) - 4 = {rate_of_change:.6f}")

# The final answer in the requested format
final_answer = rate_of_change
print(f"\n<<<{final_answer:.6f}>>>")