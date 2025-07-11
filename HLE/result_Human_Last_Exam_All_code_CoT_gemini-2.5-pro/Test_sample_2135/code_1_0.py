import numpy as np
from scipy.integrate import solve_ivp

def pendulum_motion(t, y):
  """
  Defines the differential equation for the pendulum's motion: dy/dt = sin(y(t) - t) - 4.
  y is an array, but for this problem, it's a single element y[0].
  """
  return np.sin(y[0] - t) - 4

# Define the initial condition y(0) = 0
y0 = [0]

# Define the time at which we want to calculate the rate of change
t_final = np.pi / 6

# Numerically solve the ODE to find the value of y at t = pi/6
# We solve from t=0 to t=t_final and ask for the solution only at t_final.
solution = solve_ivp(pendulum_motion, [0, t_final], y0, t_eval=[t_final])
y_at_t_final = solution.y[0, 0]

# Now, calculate the instantaneous rate of change (dy/dt) at t = pi/6
# using the value of y(pi/6) we just found.
dydt_at_t_final = pendulum_motion(t_final, [y_at_t_final])

# Print the results, showing the numbers used in the final calculation.
print(f"To find the instantaneous rate of change dy/dt at t = π/6, we first need the value of y(π/6).")
print(f"By numerically solving the differential equation from t=0 to t=π/6, we find:")
print(f"y(π/6) ≈ {y_at_t_final:.6f}")
print("\nNow we can substitute the values into the equation for the rate of change:")
print("dy/dt = sin(y(t) - t) - 4")
print(f"dy/dt = sin({y_at_t_final:.6f} - {t_final:.6f}) - 4")
argument_of_sin = y_at_t_final - t_final
print(f"dy/dt = sin({argument_of_sin:.6f}) - 4")
value_of_sin = np.sin(argument_of_sin)
print(f"dy/dt = {value_of_sin:.6f} - 4")
print(f"dy/dt ≈ {dydt_at_t_final:.6f}")

print("\nThe instantaneous rate of change at t = π/6 is approximately -4.1818.")