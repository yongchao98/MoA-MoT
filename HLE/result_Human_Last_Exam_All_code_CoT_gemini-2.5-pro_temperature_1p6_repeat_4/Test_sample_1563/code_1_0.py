import numpy as np
from scipy.integrate import solve_ivp

# Define the system of differential equations
# y[0] represents a(t), and y[1] represents b(t)
def ode_system(t, y, A, k):
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# Set the given parameters
A = 1
k = 5

# Set the initial conditions
a0 = 0.1
b0 = 2
y0 = [a0, b0]

# Define the time interval for the integration.
# We choose a sufficiently large interval to ensure the event is found.
t_span = [0, 20]

# Define the event function to find when b(t) = 0.5.
# The solver will find the root of this function.
def event_b_is_half(t, y, A, k):
    return y[1] - 0.5

# We want the integration to stop when the event occurs (terminal=True).
# Since b(t) starts at 2 and we are looking for 0.5, we expect its value
# to be decreasing, so we set the direction of the zero-crossing to -1.
event_b_is_half.terminal = True
event_b_is_half.direction = -1

# Solve the ODE system
solution = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=y0,
    args=(A, k),
    events=event_b_is_half,
    dense_output=True # Use dense output for more accurate event location
)

# Output the results, including the numbers in the final equations
print("The system of differential equations with the given parameters is:")
# The following line prints the first equation with its numbers
print(f"a'(t) = -0.5 * a^2(t) - {A} * b^2(t) + {k} * (b(t) - 1)")
# The following line prints the second equation
print(f"b'(t) = -a(t) * b(t)")
print(f"\nThe initial conditions are (a(0), b(0)) = ({a0}, {b0}).")

# The final equation we are solving for is b(t) = 0.5
target_b = 0.5
print(f"\nWe need to find the time t for the final equation: b(t) = {target_b}")

# Check if the event was found and print the time
if solution.t_events and solution.t_events[0].size > 0:
    time_of_event = solution.t_events[0][0]
    print(f"\nThe time t at which b(t) = {target_b} is approximately: {time_of_event:.4f}")
    
    # Determine which answer choice the result falls into
    if 1 < time_of_event < 2:
      print("This time is in the range (1, 2).")
    elif 2 < time_of_event < 6:
      print("This time is in the range (2, 6).")
    elif 6 < time_of_event < 10:
      print("This time is in the range (6, 10).")
    elif time_of_event > 10:
      print("This time is greater than 10.")
    else:
      print("This time is not in the specified ranges.")
else:
    print("\nThe condition b(t) = 0.5 was not met in the integration time span.")
