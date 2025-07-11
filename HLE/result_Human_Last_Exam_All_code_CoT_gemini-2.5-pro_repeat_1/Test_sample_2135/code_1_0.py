import numpy as np
from scipy.integrate import solve_ivp

# Define the differential equation for the pendulum's motion
def pendulum_motion(t, y):
    """
    Represents the differential equation dy/dt = sin(y(t) - t) - 4.
    y is expected to be an array-like object, so we access its first element y[0].
    """
    return np.sin(y[0] - t) - 4

# Set the initial condition y(0) = 0
y0 = [0]

# Set the time at which we want to find the rate of change
t_eval = np.pi / 6

# Solve the ordinary differential equation (ODE)
# We solve from t=0 to t=pi/6 to find the value of y at pi/6
# dense_output=True allows us to get a continuous solution
solution = solve_ivp(
    pendulum_motion,
    [0, t_eval],
    y0,
    dense_output=True
)

# Extract the value of y at t = pi/6
y_at_t_eval = solution.sol(t_eval)[0]

# Calculate the instantaneous rate of change (dy/dt) at t = pi/6
# by plugging the values of t and y(t) into the original equation
rate_of_change = pendulum_motion(t_eval, [y_at_t_eval])

# Print the steps and the final result
print(f"The task is to calculate dy/dt at t = π/6 for the equation dy/dt = sin(y(t) - t) - 4.")
print(f"Given the initial condition y(0) = 0.\n")
print(f"Step 1: Find the value of y(t) at t = π/6.")
print(f"Using a numerical solver, for t = π/6 ≈ {t_eval}, we get:")
print(f"y(π/6) ≈ {y_at_t_eval}\n")

print(f"Step 2: Substitute t and y(t) into the equation to find dy/dt.")
print(f"dy/dt = sin(y(π/6) - π/6) - 4")
print(f"dy/dt = sin({y_at_t_eval} - {t_eval}) - 4")
y_minus_t = y_at_t_eval - t_eval
print(f"dy/dt = sin({y_minus_t}) - 4")
sin_val = np.sin(y_minus_t)
print(f"dy/dt = {sin_val} - 4\n")

print(f"The final calculated instantaneous rate of change at t = π/6 is:")
print(rate_of_change)
