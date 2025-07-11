import numpy as np
from scipy.integrate import solve_ivp

# 1. Define the system of differential equations
def system_of_odes(t, y, A, k):
    """
    Defines the system of ODEs.
    y[0] = a(t), y[1] = b(t)
    """
    a, b = y
    # Equation for a'(t)
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    # Equation for b'(t)
    dbdt = -a * b
    return [dadt, dbdt]

# 2. Define the event function to find when b(t) = 0.5
def event_b_equals_0_5(t, y, A, k):
    """
    Event function that triggers when b(t) - 0.5 = 0.
    """
    return y[1] - 0.5

# Make the integration stop when the event is found
event_b_equals_0_5.terminal = True
# We are looking for b(t) to be decreasing, so we set the direction
event_b_equals_0_5.direction = -1

# 3. Set the parameters and initial conditions
k = 5
A = 1
a0 = 0.1
b0 = 2
initial_conditions = [a0, b0]

# 4. Set the time span for the integration
# We choose a sufficiently long time interval to ensure the event is found.
t_span = [0, 20]

print("Solving the system of differential equations:")
print("a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)")
print("b'(t) = -a*b")
print(f"With parameters: k = {k}, A = {A}")
print(f"And initial conditions: a(0) = {a0}, b(0) = {b0}")
print("\nFinding the time t when b(t) = 0.5...")

# 5. Solve the ODE system
solution = solve_ivp(
    fun=system_of_odes,
    t_span=t_span,
    y0=initial_conditions,
    args=(A, k),
    events=event_b_equals_0_5,
    dense_output=True
)

# 6. Extract and print the result
# The times of the events are stored in the 't_events' attribute
if solution.t_events[0].size > 0:
    event_time = solution.t_events[0][0]
    print(f"\nThe condition b(t) = 0.5 is met at t = {event_time:.4f}")
else:
    print("\nThe condition b(t) = 0.5 was not met in the specified time interval.")
