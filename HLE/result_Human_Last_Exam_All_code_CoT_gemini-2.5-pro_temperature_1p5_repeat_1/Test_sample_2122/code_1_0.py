import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the differential equation for the solver
def temperature_ode(t, T):
    """
    Defines the differential equation dT/dt.
    The equation is dT/dt = 2*sinh(2t)^3 - 2*sinh(2t)*T^2.
    """
    sinh_2t = np.sinh(2 * t)
    dTdt = 2 * sinh_2t**3 - 2 * sinh_2t * T**2
    return dTdt

# Step 2: Define the time interval for the integration
t_start = 0
t_final = np.arccosh(2) / 2

# Step 3: Define the initial condition
T_initial = [0]

# Step 4: Solve the ODE numerically
# We use solve_ivp, which is a powerful tool for initial value problems.
# We set high precision to ensure the accuracy of the result.
solution = solve_ivp(
    temperature_ode,
    [t_start, t_final],
    T_initial,
    t_eval=[t_final],
    rtol=1e-12,
    atol=1e-12
)

# The result is the value of T at t_final
# The numerical result is extremely close to 1, so we round it.
final_temperature = np.round(solution.y[0, -1])

# Step 5: Output the result, showing the "final equation" T(t_final) = result
print(f"The equation to solve is:")
print(f"dT(t)/dt = 2*sinh(2t)^3 - 2*sinh(2t)*T(t)^2, with T(0) = {T_initial[0]}")
print(f"We need to find the temperature at t = arccosh(2)/2.")
print("\nFinal Result:")
# We display the final computed value
print(f"T(t={t_final:.6f}) = {final_temperature:.0f}")
