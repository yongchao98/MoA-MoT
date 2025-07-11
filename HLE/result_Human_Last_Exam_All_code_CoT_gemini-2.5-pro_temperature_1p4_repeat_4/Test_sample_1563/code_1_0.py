import numpy as np
from scipy.integrate import solve_ivp

# Define the constants for the system of differential equations
k = 5
A = 1
coeff_a_sq = -0.5
coeff_b_sq = -A
coeff_b = k
const_k = -k
coeff_ab = -1

# Define the system of differential equations
# y[0] corresponds to a(t), y[1] corresponds to b(t)
def ode_system(t, y):
    a, b = y
    dadt = coeff_a_sq * a**2 + coeff_b_sq * b**2 + coeff_b * b + const_k
    dbdt = coeff_ab * a * b
    return [dadt, dbdt]

# Set the initial conditions
a0 = 0.1
b0 = 2
initial_conditions = [a0, b0]

# We need to find the time 't' when b(t) = 0.5.
# We define an event function that the solver will monitor.
# The event occurs when this function's value is zero.
def event_b_reaches_target(t, y):
    # y[1] is b(t), so we are looking for y[1] - 0.5 = 0
    return y[1] - 0.5

# We want the integration to stop when the event is found.
event_b_reaches_target.terminal = True
# Since b(t) starts at 2 and we are looking for it to reach 0.5,
# the value of b(t) is decreasing. So, the direction of the zero-crossing is -1.
event_b_reaches_target.direction = -1

# Set a sufficiently long time span for the integration to ensure the event is found.
t_span = [0, 20]

# Use solve_ivp to find the solution.
solution = solve_ivp(
    ode_system,
    t_span,
    initial_conditions,
    events=event_b_reaches_target,
    dense_output=True # Use dense_output for accurate event location
)

# Output the results
print("Solving the system of differential equations:")
print(f"a'(t) = ({coeff_a_sq})*a(t)^2 + ({coeff_b_sq})*b(t)^2 + ({coeff_b})*b(t) + ({const_k})")
print(f"b'(t) = ({coeff_ab})*a(t)*b(t)")
print(f"\nWith initial conditions: a(0) = {a0}, b(0) = {b0}")

# Check if the event occurred and extract the time.
if solution.t_events[0].size > 0:
    time_of_event = solution.t_events[0][0]
    print(f"\nThe time t at which b(t) = 0.5 is estimated to be: {time_of_event:.4f}")

    if 1 < time_of_event < 2:
        answer_choice = 'A'
    elif 2 < time_of_event < 6:
        answer_choice = 'B'
    elif 6 < time_of_event < 10:
        answer_choice = 'C'
    elif time_of_event >= 10:
        answer_choice = 'D'
    else: # Should not happen with these ranges
        answer_choice = 'E'
    print(f"This time falls into the range for answer choice {answer_choice}.")
else:
    print("\nThe condition b(t) = 0.5 was not met in the specified time interval.")
    answer_choice = 'E'
    print(f"This corresponds to answer choice {answer_choice}.")

<<<C>>>