import numpy as np
from scipy.integrate import solve_ivp

# 1. Define the system of differential equations
def system_of_des(t, y, A, k):
    """
    Defines the system of ODEs.
    y[0] represents a(t)
    y[1] represents b(t)
    """
    a, b = y
    # a'(t) = -1/2*a^2(t) - A*b^2(t) + k*(b(t)-1)
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    # b'(t) = -a(t)*b(t)
    dbdt = -a * b
    return [dadt, dbdt]

# 2. Set parameters and initial conditions
A = 1.0
k = 5.0
a0 = 0.1
b0 = 2.0
y0 = [a0, b0]

# 3. Define the event function to find when b(t) = 0.5
def event_b_is_0_5(t, y, A, k):
    """
    Event function that returns 0 when b(t) = 0.5.
    """
    return y[1] - 0.5

# We expect b(t) to decrease from 2 to 0.5, so the event happens when the function's value crosses zero from positive to negative.
event_b_is_0_5.direction = -1
# We want to stop the integration once the event is found.
event_b_is_0_5.terminal = True

# 4. Set the time span for the integration.
# We choose a sufficiently large interval to ensure the event is captured.
t_span = (0, 20)

# 5. Solve the ODE system
sol = solve_ivp(
    fun=system_of_des,
    t_span=t_span,
    y0=y0,
    args=(A, k),
    events=event_b_is_0_5,
    dense_output=True  # Recommended for accurate event location
)

# 6. Print the results
print("The system of differential equations with the given parameters is:")
# As requested, outputting the numbers in the final equation
print(f"a'(t) = -0.5 * a(t)^2 - {A} * b(t)^2 + {k} * (b(t) - 1)")
print(f"b'(t) = -a(t) * b(t)")
print(f"\nWith initial conditions a(0) = {a0} and b(0) = {b0}.")

# Check if the event was found and print the time
if sol.status == 1 and sol.t_events[0].size > 0:
    t_event = sol.t_events[0][0]
    print(f"\nThe time 't' at which b(t) = 0.5 is estimated to be:")
    print(t_event)
else:
    print("\nThe event b(t) = 0.5 was not found in the specified time interval.")
