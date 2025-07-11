import numpy as np
from scipy.integrate import solve_ivp

# 1. Define the system of ODEs with given parameters
# The system is:
# a'(t) = -0.5*a^2(t) - A*b^2(t) + k*(b(t)-1)
# b'(t) = -a(t)*b(t)
# We are given A=1 and k=5.
# y[0] represents a(t) and y[1] represents b(t).
def ode_system(t, y):
    a, b = y
    A = 1.0
    k = 5.0
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    db_dt = -a * b
    return [da_dt, db_dt]

# 2. Define the event function to find when b(t) = 0.5
# The solver will look for a root of this function.
def event_b_reaches_target(t, y):
    return y[1] - 0.5

# This option makes the integration stop when the event is found.
event_b_reaches_target.terminal = True
# This option specifies the direction of the zero crossing; b(t) is decreasing.
event_b_reaches_target.direction = -1

# 3. Set the initial conditions
a0 = 0.1
b0 = 2.0
initial_state = [a0, b0]

# 4. Set the time span for integration.
# We set a large upper limit, as the solver will stop when the event is found.
t_span = [0, 100]

# 5. Solve the ODE system
solution = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=initial_state,
    events=event_b_reaches_target,
    dense_output=True
)

# 6. Output the final result
print("Solving the system of differential equations:")
# As requested, output each number in the final equation
print("a'(t) = -0.5*a(t)^2 - 1*b(t)^2 + 5*(b(t) - 1)")
print("b'(t) = -1*a(t)*b(t)")
print(f"with initial conditions a(0) = {a0} and b(0) = {b0}.")
print("\nWe are looking for the time t when b(t) = 0.5.")

if solution.status == 1 and solution.t_events[0].size > 0:
    # solution.status == 1 means a terminal event was found.
    t_at_event = solution.t_events[0][0]
    print(f"\nThe estimated time t is: {t_at_event:.4f}")
else:
    print("\nA time t where b(t) = 0.5 was not found in the integration interval.")
