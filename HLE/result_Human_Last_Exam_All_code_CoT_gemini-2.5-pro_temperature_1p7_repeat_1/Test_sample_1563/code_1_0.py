import numpy as np
from scipy.integrate import solve_ivp

# --- Step 1: Define parameters and initial conditions ---
k = 5.0
A = 1.0
a0 = 0.1
b0 = 2.0
y0 = [a0, b0]
b_target = 0.5

# --- Step 2: Define the ODE system for the solver ---
# The state vector y is [a(t), b(t)]
def ode_system(t, y):
    """Defines the system of differential equations."""
    a, b = y
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    db_dt = -a * b
    return [da_dt, db_dt]

# --- Step 3: Define the event to find when b(t) = 0.5 ---
def find_b_target(t, y):
    """Event function to find t when b(t) equals b_target."""
    return y[1] - b_target

# Make the integration stop when the event occurs
find_b_target.terminal = True
# Trigger event when b(t) is decreasing through the target
find_b_target.direction = -1

# --- Step 4: Print problem statement and solve the ODE system ---
print("Solving the system of differential equations:")
# As requested, output each number in the final equation.
# a'(t) = -0.5*a^2 - A*b^2 + k*b - k
print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*b(t) - {k}")
print(f"b'(t) = -a(t)*b(t)")
print(f"with initial conditions a(0) = {a0}, b(0) = {b0}")
print("-" * 30)

# Set a sufficiently long time span to ensure the event is found.
t_span = (0, 20)

# Use solve_ivp to integrate the system.
solution = solve_ivp(
    ode_system,
    t_span,
    y0,
    events=find_b_target,
    dense_output=True # Use for accurate event location
)

# --- Step 5: Output the result ---
# The time of the event is stored in the 't_events' attribute.
if solution.status == 1 and len(solution.t_events[0]) > 0:
    t_at_event = solution.t_events[0][0]
    print(f"The time t at which b(t) reaches {b_target} is: {t_at_event:.4f}")
else:
    print(f"The condition b(t) = {b_target} was not met in the time interval {t_span}.")
    # Additional info can be checked from the solution object if needed
    # print(f"Solver message: {solution.message}")
    # print(f"Integration ended at t={solution.t[-1]} with b(t)={solution.y[1][-1]}")

<<<B>>>