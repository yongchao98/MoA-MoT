import numpy as np
from scipy.integrate import solve_ivp

# --- Problem Setup and Solution ---

# 1. Define the parameters and initial conditions for the system.
# a'(t) = -1/2 * a^2 - A * b^2 + k * (b - 1)
# b'(t) = -a * b
k = 5.0
A = 1.0
a0 = 0.1
b0 = 2.0
initial_conditions = [a0, b0]

# 2. Define the ODE system as a Python function.
def ode_system(t, y):
    """Defines the system of differential equations."""
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# 3. Define an event to find the time when b(t) = 0.5.
# The solver will find the root of this function.
def event_b_target(t, y):
    """Event function that triggers when b(t) = 0.5."""
    return y[1] - 0.5

# Stop the integration when the event occurs.
event_b_target.terminal = True
# Since b(t) starts at 2 and is expected to decrease to 0.5, 
# we specify the direction of the zero-crossing.
event_b_target.direction = -1

# 4. Solve the initial value problem.
# We set a sufficiently large time span to ensure the event is found.
t_span = [0, 20]

solution = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=initial_conditions,
    events=event_b_target,
    dense_output=True # Use dense output for accurate event location.
)

# 5. Print the results.
print("We are solving the following system of differential equations:")
print(f"a'(t) = -0.5*a(t)^2 - {A:.1f}*b(t)^2 + {k:.1f}*(b(t)-1)")
print(f"b'(t) = -a(t)*b(t)")
print(f"with initial conditions a(0) = {a0} and b(0) = {b0}.")
print("\nWe want to find the time t at which b(t) = 0.5.")
print("-" * 50)

# Check if the event was found and report the time.
if solution.t_events[0].size > 0:
    time_of_event = solution.t_events[0][0]
    print(f"The condition b(t) = 0.5 is met at t = {time_of_event:.4f}")
else:
    print("A time t where b(t) = 0.5 was not found in the integration interval.")
