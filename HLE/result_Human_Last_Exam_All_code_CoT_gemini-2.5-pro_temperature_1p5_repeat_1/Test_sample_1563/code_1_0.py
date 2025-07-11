import numpy as np
from scipy.integrate import solve_ivp

# 1. Define the system of ODEs with the given parameters
# a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
# b'(t) = -a*b
def ode_system(t, y, A, k):
    """
    Defines the system of differential equations.
    y[0] = a(t), y[1] = b(t)
    """
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# 2. Define the event function to find when b(t) = 0.5
def event_b_reaches_half(t, y, A, k):
    """
    Event function that returns 0 when b(t) = 0.5.
    y[1] corresponds to b(t).
    """
    return y[1] - 0.5

# Set the event to terminate the integration when it occurs.
event_b_reaches_half.terminal = True
# The event should be triggered when b(t) is decreasing.
event_b_reaches_half.direction = -1

# 3. Set parameters and initial conditions
k = 5.0
A = 1.0
a0 = 0.1
b0 = 2.0
y0 = [a0, b0]

# Define the time interval for the integration. A large enough interval is chosen.
t_span = (0, 30)

# Print the problem being solved, including all numbers in the equations.
print("Solving the system of differential equations:")
print(f"a'(t) = -0.5 * a(t)**2 - {A} * b(t)**2 + {k} * (b(t) - 1)")
print(f"b'(t) = -a(t) * b(t)")
print("\nWith initial conditions:")
print(f"a(0) = {a0}")
print(f"b(0) = {b0}")
print("\nObjective: Find the time t at which b(t) = 0.5")
print("\nSolving...")

# 4. Solve the ODE system
sol = solve_ivp(
    ode_system,
    t_span,
    y0,
    args=(A, k),
    events=event_b_reaches_half,
    dense_output=True  # For a smooth solution curve
)

# 5. Output the result
print("\nResult:")
# Check if the event was successfully found
if sol.status == 1 and len(sol.t_events[0]) > 0:
    t_event = sol.t_events[0][0]
    print(f"The time t at which b(t) = 0.5 is approximately: {t_event:.4f}")
else:
    # This case corresponds to choice E
    print("A time t where b(t) = 0.5 was not found in the integration interval.")
    print("This might happen if no such t exists or the time span is too short.")
