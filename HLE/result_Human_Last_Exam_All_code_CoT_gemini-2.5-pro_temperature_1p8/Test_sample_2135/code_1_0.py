import numpy as np
from scipy.integrate import solve_ivp

def pendulum_motion(t, y):
    """
    Defines the differential equation for the pendulum's motion:
    dy/dt = sin(y(t) - t) - 4
    """
    # The input y is an array, so we access its first element.
    return np.sin(y[0] - t) - 4

# Define the time point of interest and the initial condition.
t_final = np.pi / 6
y0 = [0]  # Initial condition y(0) = 0

# We need to solve the ODE to find the value of y at t = pi/6.
# The time span for the integration is from 0 to pi/6.
t_span = [0, t_final]

# Use a numerical solver (solve_ivp) to find the solution.
# 'dense_output=True' allows us to get a continuous solution.
solution = solve_ivp(pendulum_motion, t_span, y0, dense_output=True)

# From the solution, get the value of y at t = pi/6.
y_at_t_final = solution.sol(t_final)[0]

# Now, calculate the instantaneous rate of change (dy/dt) at t = pi/6
# using the value of y we just found.
rate_of_change = pendulum_motion(t_final, [y_at_t_final])

# Display the step-by-step calculation.
print(f"To find the instantaneous rate of change dy/dt at t = π/6, we use the equation:")
print(f"dy/dt = sin(y(t) - t) - 4")
print("\nFirst, we solve the initial value problem to find y(π/6) given y(0) = 0.")
print(f"Numerically solving the ODE from t=0 to t=π/6 gives y(π/6) ≈ {y_at_t_final:.6f}")
print("\nNext, we substitute t = π/6 and the calculated y(π/6) into the equation:")
# We use np.pi/6 for the symbolic representation and t_final for the value.
print(f"dy/dt |_(t=π/6) = sin(y(π/6) - π/6) - 4")
print(f"dy/dt = sin({y_at_t_final:.6f} - {t_final:.6f}) - 4")
argument_of_sin = y_at_t_final - t_final
print(f"dy/dt = sin({argument_of_sin:.6f}) - 4")
sin_value = np.sin(argument_of_sin)
print(f"dy/dt = {sin_value:.6f} - 4")
print(f"dy/dt = {rate_of_change:.6f}")
print("\nTherefore, the final result is:")
print(f"{rate_of_change}")