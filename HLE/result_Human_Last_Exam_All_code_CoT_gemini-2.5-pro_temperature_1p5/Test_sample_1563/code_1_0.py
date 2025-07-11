import numpy as np
from scipy.integrate import solve_ivp

# This script solves a system of ordinary differential equations (ODEs)
# to find the time t at which b(t) reaches a specific value.

# --- Problem Definition ---
# The system of ODEs is:
# a'(t) = -0.5 * a(t)^2 - A * b(t)^2 + k * (b(t) - 1)
# b'(t) = -a(t) * b(t)

# Parameters and Initial Conditions:
A = 1
k = 5
a0 = 0.1
b0 = 2
b_target = 0.5

# We print the equation with the given values.
print("Solving the system of differential equations:")
print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
print(f"b'(t) = -a(t)*b(t)")
print(f"with initial conditions a(0) = {a0}, b(0) = {b0}.")
print(f"The goal is to find the time t at which b(t) = {b_target}.")
print("-" * 30)


# --- Numerical Solution ---

# 1. Define the ODE system as a function for the solver.
# The state vector y is [a(t), b(t)].
def ode_system(t, y, A_param, k_param):
    """Defines the system of differential equations."""
    a, b = y
    dadt = -0.5 * a**2 - A_param * b**2 + k_param * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# 2. Define an event function to find when b(t) reaches the target value.
# The solver will find the root of this function.
def find_b_target(t, y, A_param, k_param):
    """Event function: returns b(t) - b_target."""
    return y[1] - b_target

# The integration should stop when this event occurs.
find_b_target.terminal = True
# We expect b(t) to be decreasing, so we look for a negative-going zero-crossing.
find_b_target.direction = -1

# 3. Set up the solver inputs.
# Initial conditions vector
y0 = [a0, b0]
# Time interval for integration. Should be long enough for the event to occur.
t_span = [0, 20]

# 4. Call the ODE solver.
# We pass the ODE system, time span, initial conditions, and the event function.
# The args parameter is used to pass the constants A and k to our functions.
solution = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=y0,
    events=find_b_target,
    args=(A, k),
    dense_output=True
)

# 5. Process and print the results.
# Check if the event was found. sol.t_events will be non-empty.
if solution.t_events[0].size > 0:
    # The time at which the event occurred
    t_event = solution.t_events[0][0]
    print(f"Successfully found the event.")
    print(f"The time t at which b(t) = {b_target} is: {t_event:.4f}")

    # Determine which answer choice interval the time falls into.
    if 1 < t_event < 2:
        answer_choice = "A"
    elif 2 < t_event < 6:
        answer_choice = "B"
    elif 6 < t_event < 10:
        answer_choice = "C"
    elif t_event > 10:
        answer_choice = "D"
    else:
        # This case is for completeness, but not expected based on options.
        answer_choice = "E"
elif solution.status == 0:
    # This means the integration finished the whole t_span without the event.
    print(f"The integration finished at t = {solution.t[-1]} without b(t) reaching {b_target}.")
    print(f"The final value of b(t) was {solution.y[1, -1]:.4f}.")
    answer_choice = "E" # No such t exists
else:
    # The solver encountered an error.
    print(f"Solver failed with status {solution.status}: {solution.message}")
    answer_choice = None # Indicate error

# Final Answer Output
if answer_choice:
    print(f"The time {t_event:.4f} is in the range 6 < t < 10.")
    print(f"<<<{answer_choice}>>>")
else:
    print("Could not determine the answer.")
