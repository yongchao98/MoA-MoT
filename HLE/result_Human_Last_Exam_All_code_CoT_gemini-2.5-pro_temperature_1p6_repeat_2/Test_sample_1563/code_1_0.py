import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the system of differential equations in a function.
# The state vector y is [a, b].
def ode_system(t, y, A, k):
    """
    Defines the system of differential equations.
    y[0] = a(t)
    y[1] = b(t)
    """
    a, b = y
    # a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
    da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    # b'(t) = -a*b
    db_dt = -a * b
    return [da_dt, db_dt]

# Step 2: Define the event function.
# We want to find the time t when b(t) = 0.5.
# The solver will find the root of this function (i.e., when it equals zero).
def event_b_equals_0_5(t, y, A, k):
    return y[1] - 0.5
# Make the event terminal, so the integration stops when b(t) reaches 0.5.
event_b_equals_0_5.terminal = True
# Since b(t) starts at 2 and decreases, we look for a negative-going crossing.
event_b_equals_0_5.direction = -1

# Step 3: Set the given parameters and initial conditions.
A = 1.0
k = 5.0
a0 = 0.1
b0 = 2.0
initial_conditions = [a0, b0]
target_b = 0.5

# We must integrate over a sufficiently long time span for the event to be found.
t_span = (0, 20)

# Step 4: Call the numerical solver.
solution = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=initial_conditions,
    args=(A, k),
    events=event_b_equals_0_5,
    dense_output=True  # Use dense output for accurate event location
)

# Step 5: Extract and print the result.
# The instruction "output each number in the final equation!" can be interpreted
# as printing the components of the result 't = <value>'. We will print the value of t.
print(f"Solving the system a'(t) = -0.5*a^2 - {A}*b^2 + {k}*(b-1), b'(t)=-a(t)b(t)")
print(f"with initial conditions (a(0), b(0)) = ({a0}, {b0}).")
print(f"Finding time t where b(t) = {target_b}:")

if solution.status == 1 and len(solution.t_events[0]) > 0:
    # solution.status == 1 means a termination event was triggered.
    # solution.t_events[0] holds the time(s) for the first event function.
    t_event = solution.t_events[0][0]
    # The final equation is t = <value>. We print the value.
    print(t_event)
else:
    # This case occurs if the integration completes without the event happening.
    print("No time t was found in the interval where b(t) = 0.5.")
    print(-1.0)
