import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the system of differential equations with the given parameters.
# The system is:
# a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
# b'(t) = -a*b
# with A=1, k=5.

# Parameters for the equation
A = 1
k = 5

def differential_equations(t, y):
    """
    Defines the system of ODEs.
    
    Args:
        t: time
        y: a list or array [a(t), b(t)]
        
    Returns:
        A list of the derivatives [a'(t), b'(t)]
    """
    a, b = y
    # The equations with the specified parameter values
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    db_dt = -a * b
    return [da_dt, db_dt]

# Step 2: Set the initial conditions and time span for the integration.
# Initial conditions: a(0) = 0.1, b(0) = 2
initial_conditions = [0.1, 2]

# Time interval for integration. We choose a sufficiently large end time
# to ensure the condition can be met.
t_span = (0, 20)

# Step 3: Define the event to find when b(t) = 0.5
def find_b_equals_0_5(t, y):
    """
    Event function that triggers when b(t) = 0.5.
    The solver will find the root of this function.
    """
    return y[1] - 0.5

# Make the integration stop when the event is found.
find_b_equals_0_5.terminal = True
# We expect b(t) to be decreasing from 2 towards 0.5.
find_b_equals_0_5.direction = -1

# Step 4: Solve the ODE system numerically.
solution = solve_ivp(
    fun=differential_equations,
    t_span=t_span,
    y0=initial_conditions,
    events=find_b_equals_0_5,
    dense_output=True # This provides a more accurate event location.
)

# Step 5: Print the results.
print("Problem Definition:")
print("System of differential equations:")
# This print statement fulfills the requirement to "output each number in the final equation"
print(f"  a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
print(f"  b'(t) = -a(t)*b(t)")
print(f"Initial conditions: a(0) = {initial_conditions[0]}, b(0) = {initial_conditions[1]}")
print("Objective: Find the time t at which b(t) = 0.5\n")

# Check if the event was found and print the time.
if solution.status == 1 and solution.t_events[0].size > 0:
    event_time = solution.t_events[0][0]
    print(f"Result:")
    print(f"The condition b(t) = 0.5 was met at t = {event_time:.4f}")
    
    if 1 < event_time < 2:
        choice = 'A'
    elif 2 < event_time < 6:
        choice = 'B'
    elif 6 < event_time < 10:
        choice = 'C'
    elif event_time > 10:
        choice = 'D'
    else:
        choice = "Outside of given ranges."
    print(f"This time falls into the range: {choice}")
elif solution.status == 0:
    print("Result:")
    print(f"The integration completed up to t = {solution.t[-1]} without b(t) reaching 0.5.")
    print("This corresponds to choice: E")
else:
    print(f"An issue occurred during integration. Solver status: {solution.message}")

<<<C>>>