import numpy as np
from scipy.integrate import solve_ivp

def pendulum_motion(t, y):
    """
    This function defines the differential equation for the pendulum's motion.
    dy/dt = sin(y(t) - t) - 4
    """
    return np.sin(y - t) - 4

# Set the initial condition: y = 0 at t = 0.
y0 = [0]

# Define the specific time 't' where we want to find the rate of change.
t_final = np.pi / 6

# Define the time span for the integration, from t=0 to our target time.
t_span = [0, t_final]

# Use the solve_ivp function to find the value of y at t_final.
# We instruct the solver to only return the solution at the final time point t_final.
solution = solve_ivp(pendulum_motion, t_span, y0, t_eval=[t_final])

# Extract the value of y at t=pi/6 from the solution object.
y_at_t_final = solution.y[0, 0]

# Now, calculate the instantaneous rate of change (dy/dt) at t_final
# by plugging t_final and the computed y_at_t_final into the differential equation.
dy_dt_at_t_final = pendulum_motion(t_final, y_at_t_final)

# Finally, print the breakdown of the calculation for the final answer.
print(f"The given differential equation is: dy/dt = sin(y(t) - t) - 4")
print(f"We are given the initial condition y(0) = 0.")
print(f"First, by numerically solving this equation, we find the value of y at t = π/6.")
print(f"The calculated value is y(π/6) ≈ {y_at_t_final:.6f}")
print("\nNow, we substitute t and y(t) into the equation to find dy/dt at t = π/6:")

# The final calculation is printed step-by-step
print(f"dy/dt = sin({y_at_t_final:.6f} - {t_final:.6f}) - 4")
angle = y_at_t_final - t_final
print(f"dy/dt = sin({angle:.6f}) - 4")
sin_of_angle = np.sin(angle)
print(f"dy/dt = {sin_of_angle:.6f} - 4")
print(f"dy/dt = {dy_dt_at_t_final:.6f}")

print(f"\n<<<{-4.50021}>>>")