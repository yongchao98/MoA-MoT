import numpy as np
from scipy.integrate import solve_ivp

# Define the parameters for the equations
A = 1
k = 5

# Define the system of differential equations based on the problem statement
# a'(t) = -0.5*a(t)^2 - A*b(t)^2 + k*(b(t)-1)
# b'(t) = -a(t)*b(t)
# We will print the system with the specific values for A and k.
print("The system of differential equations with the given parameters is:")
print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
print(f"b'(t) = -a(t)*b(t)")
print("-" * 30)

def system_of_odes(t, y, A, k):
    """
    Defines the system of differential equations.
    y = [a, b]
    """
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# Set the initial conditions
a0 = 0.1
b0 = 2
y0 = [a0, b0]
print(f"Initial conditions: a(0) = {a0}, b(0) = {b0}")

# Define an event function to find the time when b(t) = 0.5
def event_b_reaches_0_5(t, y, A, k):
    """
    Event function that triggers when b(t) equals 0.5.
    The solver will find the root of this function (i.e., when it's zero).
    """
    return y[1] - 0.5

# We want the integration to stop when the event occurs.
event_b_reaches_0_5.terminal = True
# We are looking for the point where b(t) crosses 0.5 from above (decreasing).
event_b_reaches_0_5.direction = -1

# Set the time span for integration. We choose a range that we expect
# to contain the solution, e.g., [0, 20].
t_span = [0, 20]

# Use solve_ivp to find the solution and the time of the event.
solution = solve_ivp(
    system_of_odes,
    t_span,
    y0,
    args=(A, k),
    events=event_b_reaches_0_5,
    dense_output=True  # Allows for a dense solution if needed
)

# Check the solution and print the result.
if solution.status == 1 and len(solution.t_events[0]) > 0:
    # sol.status == 1 means the integration was terminated by an event.
    # sol.t_events[0] contains the time(s) of the event. We take the first one.
    time_of_event = solution.t_events[0][0]
    print(f"\nThe simulation shows that b(t) reaches 0.5 at time t = {time_of_event:.4f}.")
else:
    # This case would occur if b(t) never reaches 0.5 within the t_span.
    print("\nThe event where b(t) = 0.5 was not found in the integration interval.")
    print(f"Integration ended at t = {solution.t[-1]} with a final value of b(t) = {solution.y[1][-1]:.4f}.")
