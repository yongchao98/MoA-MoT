import numpy as np
from scipy.integrate import solve_ivp

# Define the differential equation dy/dt = f(t, y)
def pendulum_motion(t, y):
  """
  Defines the differential equation for the pendulum's motion.
  dy/dt = sin(y(t) - t) - 4
  """
  return np.sin(y - t) - 4

# Set the initial condition y(0) = 0
y0 = [0]

# Set the time for which we want to calculate the rate of change
t_final = np.pi / 6

# We need to solve the ODE from t=0 to t=t_final to find y(t_final)
t_span = [0, t_final]

# Use solve_ivp to find the solution. We ask for the solution only at t_final.
# 'dense_output=True' could also be used for a continuous solution.
sol = solve_ivp(pendulum_motion, t_span, y0, t_eval=[t_final])

# Extract the value of y at t = pi/6 from the solution
# sol.y is an array of solutions, we need the first state variable at the first evaluation time.
y_at_t_final = sol.y[0, 0]

# Now, calculate the instantaneous rate of change dy/dt at t = pi/6
# by substituting t and y(t) into the original equation.
rate_of_change = pendulum_motion(t_final, y_at_t_final)

# Print the breakdown of the calculation as requested
print("To find the instantaneous rate of change at t = π/6, we follow these steps:")
print("1. The differential equation is: dy/dt = sin(y(t) - t) - 4")
print(f"2. We first need the value of y(t) at t = π/6 ≈ {t_final:.6f}.")
print(f"3. By numerically solving the ODE with y(0) = 0, we find y(π/6) ≈ {y_at_t_final:.6f}.")
print("4. Now we substitute these values into the equation:")
print(f"   dy/dt = sin({y_at_t_final:.6f} - {t_final:.6f}) - 4")
print(f"   dy/dt = sin({y_at_t_final - t_final:.6f}) - 4")
print(f"   dy/dt = {np.sin(y_at_t_final - t_final):.6f} - 4")
print(f"5. The final calculated instantaneous rate of change is:")
print(rate_of_change)