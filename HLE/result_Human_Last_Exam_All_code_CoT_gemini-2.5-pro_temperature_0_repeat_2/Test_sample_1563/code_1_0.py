import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the system of differential equations.
# The state vector y is [a, b].
def ode_system(t, y, A, k):
    """
    Defines the system of ODEs.
    y[0] = a(t)
    y[1] = b(t)
    """
    a, b = y
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    db_dt = -a * b
    return [da_dt, db_dt]

# Step 2: Define the event function to find when b(t) = 0.5.
def event_b_reaches_0_5(t, y, A, k):
    """
    Event function that returns zero when b(t) = 0.5.
    """
    return y[1] - 0.5

# This attribute tells the solver to stop the integration when the event is found.
event_b_reaches_0_5.terminal = True
# We expect b(t) to decrease from 2 to 0.5, so we look for a negative-going zero crossing.
event_b_reaches_0_5.direction = -1

# Step 3: Set the parameters and initial conditions.
k = 5
A = 1
a0 = 0.1
b0 = 2
initial_conditions = [a0, b0]

# Define the time span for the integration.
# We choose a large enough interval to ensure the event is found.
t_span = [0, 20]

# Step 4: Solve the ODE system using solve_ivp.
# We use a lambda function to pass the parameters A and k to the ODE and event functions.
solution = solve_ivp(
    fun=lambda t, y: ode_system(t, y, A, k),
    t_span=t_span,
    y0=initial_conditions,
    events=lambda t, y: event_b_reaches_0_5(t, y, A, k),
    dense_output=True  # Recommended for accurate event location
)

# Step 5: Extract and print the result.
print("Solving the system of differential equations:")
# Print the final equations with the numerical values of the parameters.
print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
print(f"b'(t) = -a(t)*b(t)")
print(f"\nWith initial conditions a(0) = {a0} and b(0) = {b0}.")
print("We want to find the time t at which b(t) = 0.5.")

# Check if the event was found.
if solution.t_events[0].size > 0:
    # The time of the first event is stored here.
    t_at_event = solution.t_events[0][0]
    print(f"\nUsing numerical integration, the time t at which b(t) reaches 0.5 is approximately: {t_at_event:.4f}")
else:
    print("\nThe condition b(t) = 0.5 was not met within the specified time interval.")
