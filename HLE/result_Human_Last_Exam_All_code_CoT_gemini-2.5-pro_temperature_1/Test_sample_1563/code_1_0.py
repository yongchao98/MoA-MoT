import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the system of differential equations
# a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
# b'(t) = -a*b
def model(t, y, A, k):
    """
    Defines the system of ODEs.
    y[0] = a(t), y[1] = b(t)
    """
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# Step 2: Set parameters and initial conditions
k = 5
A = 1
a0 = 0.1
b0 = 2
y0 = [a0, b0]

# Step 3: Define the event to find when b(t) = 0.5
def find_b_target(t, y, A, k):
    """
    Event function that is zero when b(t) = 0.5.
    y[1] corresponds to b(t).
    """
    return y[1] - 0.5

# Stop the integration when the event is found
find_b_target.terminal = True
# The event occurs as b(t) is decreasing (from 2 to 0.5), so direction is negative.
find_b_target.direction = -1

# Step 4: Solve the ODE system numerically
# We set a time span that is long enough for the event to occur.
# Based on answer choices, a max time of 20 is safe.
t_span = (0, 20)

# Use solve_ivp to integrate the system until the event is found.
# The `args` parameter passes A and k to the model and event functions.
sol = solve_ivp(
    model,
    t_span,
    y0,
    args=(A, k),
    events=find_b_target,
    dense_output=True # dense_output is recommended for accurate event location
)

# Step 5: Output the result
print(f"Solving the system of ODEs with parameters A={A}, k={k} and initial conditions a(0)={a0}, b(0)={b0}.")

# The times of the events are stored in sol.t_events.
# We check if the event list is not empty (i.e., the event was found).
if sol.t_events[0].size > 0:
    # The first (and only, since it's terminal) event time is our answer.
    t_event = sol.t_events[0][0]
    print(f"The time t at which b(t) reaches 0.5 is approximately: {t_event:.4f}")
else:
    # This case handles if b(t) never reaches 0.5 within the specified t_span.
    print(f"The condition b(t) = 0.5 was not met in the time interval {t_span}.")

<<<C>>>