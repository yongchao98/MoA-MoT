import numpy as np
from scipy.integrate import solve_ivp

# This script solves the given system of differential equations to find the time 't'
# at which b(t) = 0.5.

# 1. Define the system of ODEs:
# a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
# b'(t) = -a*b
def ode_system(t, y, A, k):
    """
    Defines the system of ODEs.
    y[0] = a(t), y[1] = b(t)
    """
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# 2. Define the event function to find when b(t) equals 0.5
def event_b_reaches_target(t, y):
    """
    Event function that triggers when b(t) - 0.5 = 0.
    y[1] is b(t).
    """
    return y[1] - 0.5

# Configure the event: stop integration when it occurs, and detect it as b(t) decreases.
event_b_reaches_target.terminal = True
event_b_reaches_target.direction = -1

# 3. Set constants and initial conditions
k = 5.0
A = 1.0
initial_conditions = [0.1, 2.0] # [a(0), b(0)]

# 4. Set the time span for integration. A sufficiently large end time is chosen.
t_span = [0, 20]

# 5. Solve the initial value problem
solution = solve_ivp(
    fun=lambda t, y: ode_system(t, y, A, k),
    t_span=t_span,
    y0=initial_conditions,
    events=event_b_reaches_target,
    dense_output=True
)

# 6. Output the result. The "final equation" is b(t) = 0.5.
# We output the numbers involved: the time 't' we solved for and the target value '0.5'.
if solution.t_events[0].size > 0:
    time_t = solution.t_events[0][0]
    target_b = 0.5
    print(f"The condition b(t) = {target_b} is satisfied at the following time.")
    print(f"Time t: {time_t}")
    print(f"Value of b(t): {target_b}")
else:
    print("The condition b(t) = 0.5 was not met in the specified time interval.")
