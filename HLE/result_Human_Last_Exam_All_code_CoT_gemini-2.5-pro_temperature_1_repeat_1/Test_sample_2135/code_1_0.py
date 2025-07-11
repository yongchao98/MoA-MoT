import numpy as np
from scipy.integrate import solve_ivp

def pendulum_motion(t, y):
    """
    Defines the differential equation for the pendulum's motion:
    dy/dt = sin(y - t) - 4
    """
    return np.sin(y - t) - 4

# Initial condition: y(0) = 0
t_initial = 0
y_initial = [0]  # solve_ivp expects an array-like initial condition

# The time at which we want to find the rate of change
t_target = np.pi / 6

# Step 1: Solve the initial value problem (IVP) to find y(t_target).
# We integrate from t_initial to t_target.
solution = solve_ivp(
    fun=pendulum_motion,
    t_span=[t_initial, t_target],
    y0=y_initial,
    t_eval=[t_target]  # Ensure the solution is evaluated exactly at t_target
)

# Extract the value of y at t_target from the solution object.
y_at_target = solution.y[0, 0]

# Step 2: Calculate the instantaneous rate of change (dy/dt) at t_target.
# Now that we have y(t_target), we can use the original equation for dy/dt.
rate_of_change = pendulum_motion(t_target, y_at_target)

# Step 3: Print the results, showing each number in the final equation.
print(f"To find the rate of change dy/dt at t = π/6, we use the equation:")
print("dy/dt = sin(y(π/6) - π/6) - 4\n")

print(f"First, we solve for y(π/6) numerically, starting from y(0) = 0.")
print(f"The numerical value for y at t = π/6 is: y(π/6) ≈ {y_at_target:.6f}")
print(f"The value of t is: π/6 ≈ {t_target:.6f}\n")

print("Now, we substitute these numbers into the rate of change equation:")
# The prompt requires outputting each number in the final equation.
print(f"dy/dt ≈ sin({y_at_target:.6f} - {t_target:.6f}) - 4")
print(f"dy/dt ≈ sin({y_at_target - t_target:.6f}) - 4")
print(f"dy/dt ≈ {np.sin(y_at_target - t_target):.6f} - 4")
print(f"\nThe final calculated instantaneous rate of change is: {rate_of_change:.6f}")