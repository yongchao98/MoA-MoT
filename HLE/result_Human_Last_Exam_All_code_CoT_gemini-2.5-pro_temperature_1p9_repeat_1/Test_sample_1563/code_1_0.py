import numpy as np
from scipy.integrate import solve_ivp

def system_of_odes(t, y):
    """
    Defines the system of differential equations.
    y[0] = a(t)
    y[1] = b(t)
    """
    a, b = y
    
    # Constants for the system
    k = 5.0
    A = 1.0
    
    # The differential equations
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    
    return [dadt, dbdt]

def event_b_equals_half(t, y):
    """
    Event function to find when b(t) = 0.5.
    The solver will find the root of this function.
    """
    return y[1] - 0.5

# Make the integration stop when the event is found.
event_b_equals_half.terminal = True
# We are looking for the point where the function b(t) crosses 0.5 from above,
# so the direction of the event is decreasing.
event_b_equals_half.direction = -1

# Initial conditions
a0 = 0.1
b0 = 2.0
initial_conditions = [a0, b0]

# Time interval for integration. We choose a large enough end time.
t_span = (0, 30)

# Solve the ODE system
solution = solve_ivp(
    fun=system_of_odes,
    t_span=t_span,
    y0=initial_conditions,
    events=event_b_equals_half,
    dense_output=True # Ensures the event is found with high precision
)

# Check for the event and print the result
if solution.t_events[0].size > 0:
    time_of_event = solution.t_events[0][0]
    # As requested, printing the numbers in the final equation: b(t) = 0.5
    print(f"The time `t` at which b(t) = 0.5 was found.")
    print(f"Final equation: b({time_of_event:.4f}) = 0.5")
else:
    print("The condition b(t) = 0.5 was not reached in the integration time interval.")
