import numpy as np
from scipy.integrate import solve_ivp

# Define the system of differential equations as a function
# y is a vector [a, b]
def system_odes(t, y, A, k):
    """
    Defines the system of ODEs.
    a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
    b'(t) = -a*b
    """
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# Define the event function to find when b(t) = 0.5
def find_b_target(t, y, A, k):
    """
    Event function which is zero when b(t) = 0.5.
    The solver will find the root of this function.
    """
    return y[1] - 0.5

# We know b(t) starts at 2 and b'(0) is negative, so b(t) is decreasing.
# We are looking for the event when the value crosses zero from a positive value.
find_b_target.direction = -1
# We want the integration to terminate when the event is found.
find_b_target.terminal = True

# System parameters and initial conditions
A = 1
k = 5
initial_conditions = [0.1, 2] # [a(0), b(0)]

# Time span for the integration. We choose a large enough interval to ensure the event is found.
t_span = [0, 20]

# Use solve_ivp to find the solution
solution = solve_ivp(
    fun=system_odes,
    t_span=t_span,
    y0=initial_conditions,
    args=(A, k),
    events=find_b_target,
    dense_output=True # Necessary for accurate event location
)

# Extract the time of the event from the solution object
if solution.t_events[0].size > 0:
    # The event was found, so we get the time.
    time_of_event = solution.t_events[0][0]
    # To satisfy the instruction "output each number in the final equation",
    # we will print the equation t = [value]
    print(f"The time t at which b(t)=0.5 is found to be:")
    print(f"t = {time_of_event:.4f}")
else:
    # The event was not found in the specified time span
    print("The condition b(t)=0.5 was not met in the integration time span.")
