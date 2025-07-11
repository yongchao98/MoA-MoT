import numpy as np
from scipy.integrate import solve_ivp

def pendulum_ode(t, y):
    """
    Defines the differential equation for the pendulum's motion:
    dy/dt = sin(y(t) - t) - 4
    """
    return np.sin(y[0] - t) - 4

# The time point at which we want to find the rate of change
t_target = np.pi / 6

# Initial condition: y(0) = 0
y0 = [0]

# Time span for the integration from t=0 to t=t_target
t_span = [0, t_target]

# Use solve_ivp to find the value of y at t_target
# We ask the solver to return the solution specifically at our target time.
solution = solve_ivp(pendulum_ode, t_span, y0, t_eval=[t_target])

# Extract the value of y(t_target) from the solution
y_at_target = solution.y[0, 0]

# Now, calculate the instantaneous rate of change dy/dt at t_target
# using the original differential equation.
rate_of_change = pendulum_ode(t_target, [y_at_target])

# Print the step-by-step calculation
print("To find the rate of change dy/dt at t = π/6, we use the formula:")
print("dy/dt = sin(y(t) - t) - 4\n")
print(f"First, we solve for y(t) at t = π/6 ≈ {t_target:.6f}")
print(f"The numerical solution gives y(π/6) ≈ {y_at_target:.6f}\n")
print("Now we plug these values into the formula:")
print(f"dy/dt = sin({y_at_target:.6f} - {t_target:.6f}) - 4")
print(f"dy/dt = sin({y_at_target - t_target:.6f}) - 4")
print(f"dy/dt = {np.sin(y_at_target - t_target):.6f} - 4")
print(f"dy/dt = {rate_of_change:.6f}")

print("\n<<<" + str(rate_of_change) + ">>>")