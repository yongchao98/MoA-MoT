import numpy as np
from scipy.integrate import solve_ivp

def system_of_odes(t, y):
    """
    Defines the system of differential equations.
    y[0] corresponds to a(t)
    y[1] corresponds to b(t)
    """
    a, b = y
    
    # Constants from the problem description
    k = 5.0
    A = 1.0
    
    # The system of differential equations
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    db_dt = -a * b
    
    return [da_dt, db_dt]

# Initial conditions
a0 = 0.1
b0 = 2.0
initial_conditions = [a0, b0]

# Time span for the integration. We choose an interval large enough
# for the solution to evolve. The solver will stop early if the event occurs.
t_span = [0, 20]

def event_b_reaches_0_5(t, y):
    """
    Event function that triggers when b(t) = 0.5.
    The solver looks for a zero crossing of this function.
    """
    return y[1] - 0.5

# Configure the event to be terminal, i.e., stop integration when it occurs.
event_b_reaches_0_5.terminal = True
# The direction of the crossing does not matter here, but b(t) is decreasing.
event_b_reaches_0_5.direction = -1

# Solve the ODE system
solution = solve_ivp(
    fun=system_of_odes,
    t_span=t_span,
    y0=initial_conditions,
    events=event_b_reaches_0_5,
    dense_output=True  # Allows for accurate interpolation if needed
)

# Check if the event was triggered and print the result
if solution.t_events[0].size > 0:
    t_event = solution.t_events[0][0]
    print(f"The initial condition is (a(0), b(0)) = ({a0}, {b0})")
    print(f"The condition b(t) = 0.5 is met at time t = {t_event}")
else:
    print("The condition b(t) = 0.5 was not met in the specified time interval.")
