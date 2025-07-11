import numpy as np
from scipy.integrate import solve_ivp

# 1. Define the system of differential equations with the given constants.
# The system is:
# a'(t) = -0.5 * a(t)^2 - 1 * b(t)^2 + 5 * (b(t) - 1)
# b'(t) = -a(t) * b(t)

def system_of_odes(t, y):
    """
    Defines the system of ODEs.
    y[0] represents a(t), y[1] represents b(t).
    """
    a, b = y
    
    # Constants
    k = 5.0
    A = 1.0
    
    # The equations for the derivatives
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    
    return [dadt, dbdt]

# 2. Define an event function to find when b(t) = 0.5.
# The solver will find the root of this function.
def event_b_reaches_0_5(t, y):
    """Event triggers when b(t) - 0.5 = 0."""
    return y[1] - 0.5

# The integration should stop when the event is found.
event_b_reaches_0_5.terminal = True
# The event should only trigger when b(t) is decreasing.
event_b_reaches_0_5.direction = -1

# 3. Set the initial conditions and time span.
a0 = 0.1
b0 = 2.0
initial_conditions = [a0, b0]

# Set a time span that is likely long enough for the event to occur.
t_span = [0, 20]

# Print the final equation being solved
print("Solving the system of ODEs:")
print("a'(t) = -0.5*a^2 - 1*b^2 + 5*(b-1)")
print("b'(t) = -a*b")
print(f"with initial conditions a(0) = {a0}, b(0) = {b0}\n")


# 4. Solve the ODE system.
solution = solve_ivp(
    fun=system_of_odes,
    t_span=t_span,
    y0=initial_conditions,
    events=event_b_reaches_0_5,
    dense_output=True 
)

# 5. Check the results and print the time of the event.
if solution.t_events[0].size > 0:
    # The time at which b(t) reached 0.5
    t_event = solution.t_events[0][0]
    print(f"The condition b(t) = 0.5 is met at time t = {t_event:.4f}")

    # Determine which answer choice is correct
    if 1 < t_event < 2:
        print("This time falls in the range 1 < t < 2.")
    elif 2 < t_event < 6:
        print("This time falls in the range 2 < t < 6.")
    elif 6 < t_event < 10:
        print("This time falls in the range 6 < t < 10.")
    elif t_event > 10:
        print("This time falls in the range t > 10.")
else:
    print(f"The event b(t) = 0.5 was not found in the time interval up to t = {t_span[1]}.")
    print("This implies that no such t exists or a longer time interval is needed.")

<<<C>>>