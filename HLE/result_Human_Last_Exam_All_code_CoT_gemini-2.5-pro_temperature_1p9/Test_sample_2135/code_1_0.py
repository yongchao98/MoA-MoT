import numpy as np
from scipy.integrate import solve_ivp

# Define the differential equation for the pendulum's motion
def pendulum_ode(t, y):
    """
    This function represents the differential equation dy/dt = sin(y(t) - t) - 4.
    It takes time t and the angle y as input and returns the rate of change dy/dt.
    """
    return np.sin(y - t) - 4

# Define the parameters for the problem
# Initial condition: y = 0 at t = 0
initial_condition = [0]
# The specific time at which we want to find the rate of change
t_target = np.pi / 6
# The time span for numerical integration is from t=0 to t=pi/6
t_span = [0, t_target]

print("--- Step 1: Solving for the angle y(t) at t = pi/6 ---")
# Use the `solve_ivp` function from the SciPy library to solve the initial value problem.
# We request a 'dense' output to get a continuous solution, which we can evaluate at any point.
solution = solve_ivp(
    pendulum_ode,
    t_span,
    initial_condition,
    dense_output=True,
    method='RK45'
)

# Evaluate the solution function at our target time to find y(pi/6)
y_at_target = solution.sol(t_target)[0]

print(f"The value of t is pi/6, which is approximately {t_target:.6f} seconds.")
print(f"The calculated value of the angle y at this time is y(pi/6) ≈ {y_at_target:.6f} radians.")

print("\n--- Step 2: Calculating the instantaneous rate of change dy/dt ---")
# Now, calculate the rate of change using the found value of y(pi/6)
rate_of_change = pendulum_ode(t_target, y_at_target)

# Print the final calculation in a step-by-step format
print("We substitute t and y(t) into the equation: dy/dt = sin(y(t) - t) - 4")
print(f"dy/dt = sin({y_at_target:.6f} - {t_target:.6f}) - 4")

argument_of_sin = y_at_target - t_target
print(f"dy/dt = sin({argument_of_sin:.6f}) - 4")

value_of_sin = np.sin(argument_of_sin)
print(f"dy/dt = {value_of_sin:.6f} - 4")

print("\n--- Final Result ---")
print("The final calculated instantaneous rate of change is:")
print(f"{rate_of_change:.6f} rad/s")