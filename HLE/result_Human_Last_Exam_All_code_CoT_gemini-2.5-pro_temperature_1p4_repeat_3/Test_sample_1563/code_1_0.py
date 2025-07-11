import numpy as np
from scipy.integrate import solve_ivp

# Define the system of differential equations.
# The state vector y is [a(t), b(t)].
def system_of_odes(t, y):
    """
    Defines the system of ODEs:
    a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
    b'(t) = -a*b
    """
    a, b = y
    
    # Constants given in the problem
    A = 1
    k = 5
    
    # Calculate the derivatives
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    db_dt = -a * b
    
    return [da_dt, db_dt]

# Define an event function to find the time t when b(t) = 0.5.
# The solver looks for when this function returns zero.
def event_b_reaches_0_5(t, y):
    """Event function that is zero when b(t) = 0.5."""
    return y[1] - 0.5

# Configure the event to stop the integration once found.
event_b_reaches_0_5.terminal = True
# We expect b(t) to be decreasing from 2 towards 0.5, so we look for a
# negative direction crossing.
event_b_reaches_0_5.direction = -1

# Set the initial conditions
a0 = 0.1
b0 = 2
initial_conditions = [a0, b0]

# Set the time span for the integration.
# We choose a t_max of 20, which is likely sufficient to find the event.
t_span = (0, 20)

# Use solve_ivp to find the solution and the time of the event.
solution = solve_ivp(
    system_of_odes, 
    t_span, 
    initial_conditions, 
    events=event_b_reaches_0_5,
    dense_output=True # Allows for a continuous solution
)

# Output the results
print("The system of differential equations being solved is:")
print(f"a'(t) = -0.5*a(t)^2 - 1*b(t)^2 + 5*(b(t)-1)")
print(f"b'(t) = -a(t)*b(t)")
print(f"\nWith initial conditions a(0) = {a0} and b(0) = {b0}.")

# Check if the event was found and print the time.
if solution.t_events and len(solution.t_events[0]) > 0:
    t_event = solution.t_events[0][0]
    print(f"\nThe time t at which b(t) = 0.5 is: {t_event}")
else:
    print("\nThe condition b(t) = 0.5 was not met in the integration interval.")

<<<B>>>