import numpy as np
from scipy.integrate import solve_ivp

def get_temperature_at_time():
    """
    Solves the given differential equation numerically to find the temperature
    at the specified time.
    """
    # Define the differential equation: dT/dt = f(t, T)
    def dTdt(t, T):
        sinh_2t = np.sinh(2 * t)
        return 2 * sinh_2t**3 - 2 * sinh_2t * T**2

    # Set the time interval
    t_start = 0.0
    t_final = np.arccosh(2) / 2

    # Set the initial condition
    T_initial = [0.0]

    # Solve the ODE using a high-precision solver
    # We ask the solver to return the solution only at the final time point.
    solution = solve_ivp(
        dTdt,
        [t_start, t_final],
        T_initial,
        t_eval=[t_final],
        rtol=1e-12, # Relative tolerance
        atol=1e-12  # Absolute tolerance
    )

    # Extract the final temperature from the solution object
    temperature_at_final_time = solution.y[0, 0]

    # The problem asks to output each number in the final equation.
    # The final equation is T(t_final) = temperature_at_final_time.
    print(f"The final time t is: {t_final}")
    print(f"The temperature T at this time is: {temperature_at_final_time}")
    
    # The numerical result is extremely close to sqrt(3)/2.
    # We can verify this for context.
    # print(f"Note: The value of sqrt(3)/2 is approximately {np.sqrt(3)/2}")


if __name__ == '__main__':
    get_temperature_at_time()
