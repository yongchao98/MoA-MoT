import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the constants and the system of ODEs
k = 5.0
A = 1.0

def ode_system(t, y):
    """
    Defines the system of differential equations.
    y[0] represents a(t)
    y[1] represents b(t)
    """
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# Step 2: Set the initial conditions
# a(0) = 0.1, b(0) = 2
initial_conditions = [0.1, 2.0]

# Step 3: Define the event to find when b(t) = 0.5
# We are looking for the root of the event function, i.e., when event_function = 0.
def find_b_equals_0_5(t, y):
    """Event function that triggers when b(t) = 0.5"""
    return y[1] - 0.5

# Configure the event to stop the integration when it occurs.
find_b_equals_0_5.terminal = True
# The event should trigger when b(t) is decreasing towards 0.5,
# so the event function value is going from positive to negative.
find_b_equals_0_5.direction = -1

# Step 4: Solve the ODE system
# We set a large enough time span to ensure the event is found.
t_span = [0, 100]

solution = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=initial_conditions,
    events=find_b_equals_0_5,
    dense_output=True # Use dense output for accurate event location
)

# Step 5: Output the result
# The problem asks to find the time 't' that solves the equation b(t) = 0.5
target_b_value = 0.5

print(f"Solving for t in the equation: b(t) = {target_b_value}")

# Check if the event was found and print the time
if solution.t_events[0].size > 0:
    time_of_event = solution.t_events[0][0]
    print(f"The time t at which b(t) reaches {target_b_value} is approximately: {time_of_event:.4f}")
else:
    print(f"The condition b(t) = {target_b_value} was not met in the integration time span.")
