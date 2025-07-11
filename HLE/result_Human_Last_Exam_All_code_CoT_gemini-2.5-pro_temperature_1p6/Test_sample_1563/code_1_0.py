import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the constants for the problem.
A = 1.0
k = 5.0

# Step 2: Define the system of differential equations.
# The state vector y is [a, b].
def ode_system(t, y):
    a, b = y
    # a'(t) = -1/2*a^2(t) - A*b^2(t) + k*(b(t)-1)
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    # b'(t) = -a(t)*b(t)
    db_dt = -a * b
    return [da_dt, db_dt]

# Step 3: Set the initial conditions.
a0 = 0.1
b0 = 2.0
initial_conditions = [a0, b0]

# Step 4: Define an event function to find when b(t) = 0.5.
# The solver will find the time 't' where this function returns 0.
def event_b_equals_half(t, y):
    # y[1] corresponds to b(t)
    return y[1] - 0.5

# Configure the event to stop the integration when it occurs.
event_b_equals_half.terminal = True
# We are looking for b(t) to decrease to 0.5, so the direction is downward.
event_b_equals_half.direction = -1

# Step 5: Define the time span for the integration.
# We choose a span long enough to ensure the event is found.
t_span = (0, 20)

# Step 6: Solve the ODE system.
solution = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=initial_conditions,
    events=event_b_equals_half,
    dense_output=True
)

# Step 7: Output the results.
print("The system of differential equations to solve is:")
print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
print(f"b'(t) = -a(t)*b(t)")
print(f"With initial conditions a(0) = {a0} and b(0) = {b0}.")
print("\nThe final equation we need to satisfy is: b(t) = 0.5")

# Check if the event was found and print the time.
if solution.status == 1 and solution.t_events[0].size > 0:
    event_time = solution.t_events[0][0]
    # Get the value of b at the event time to confirm it's 0.5
    b_at_event_time = solution.sol(event_time)[1]
    
    print("\nResult:")
    print(f"The condition b(t) = {b_at_event_time:.1f} is satisfied at time t = {event_time:.4f}.")
else:
    print("\nThe condition b(t) = 0.5 was not met in the specified time interval.")

<<<A>>>