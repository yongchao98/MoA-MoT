import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the system of differential equations
# The state vector y is [a(t), b(t)]
def ode_system(t, y, A, k):
    """
    Defines the system of differential equations.
    a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
    b'(t) = -a*b
    """
    a, b = y
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    db_dt = -a * b
    return [da_dt, db_dt]

# Step 2: Define the event function to find when b(t) = 0.5
def event_b_reaches_0_5(t, y, A, k):
    """
    Event function that triggers when b(t) = 0.5.
    The solver will find the root of this function.
    """
    return y[1] - 0.5

# The event should stop the integration once it's found
event_b_reaches_0_5.terminal = True
# The event occurs as b(t) decreases through 0.5
event_b_reaches_0_5.direction = -1

# Step 3: Set parameters and initial conditions
A = 1.0
k = 5.0
a0 = 0.1
b0 = 2.0
y0 = [a0, b0]

# Step 4: Set the time span for integration
# We choose a time span long enough for the event to occur.
t_span = [0, 20]

# Step 5: Solve the ODE system
sol = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=y0,
    args=(A, k),
    events=event_b_reaches_0_5,
    dense_output=True # Allows for a continuous solution
)

# Step 6: Output the results
print("Solving the system of differential equations:")
print(f"a'(t) = -0.5 * a(t)^2 - {A} * b(t)^2 + {k} * (b(t) - 1)")
print(f"b'(t) = -a(t) * b(t)")
print(f"With initial conditions: a(0) = {a0}, b(0) = {b0}")
print("\nFinding the time t when b(t) = 0.5...")

# Check if the event was detected
if sol.t_events[0].size > 0:
    t_at_event = sol.t_events[0][0]
    print(f"\nThe time t at which b(t) reaches 0.5 is: {t_at_event:.4f}")
else:
    print("\nThe condition b(t) = 0.5 was not met in the integration interval.")

<<<B>>>