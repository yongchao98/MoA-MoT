import numpy as np
from scipy.integrate import solve_ivp

# 1. Define the constants and initial conditions from the problem.
k = 5
A = 1
a0 = 0.1
b0 = 2
y0 = [a0, b0]

# The target is to find the time t when b(t) equals this value.
b_target = 0.5

# 2. Define the system of differential equations.
# The function takes time t and the state vector y = [a(t), b(t)] as input.
def ode_system(t, y):
    """
    Defines the system of ODEs:
    a'(t) = -1/2*a^2 - A*b^2 + k*(b-1)
    b'(t) = -a*b
    """
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# 3. Define an event function to find when b(t) reaches the target value.
def reach_target_b(t, y):
    """
    Event function that returns 0 when b(t) = b_target.
    """
    return y[1] - b_target

# Configure the event to stop the integration when it occurs.
reach_target_b.terminal = True
# We expect b(t) to decrease from 2 to 0.5, so we are looking for a
# root where the function is decreasing (direction = -1).
reach_target_b.direction = -1

# 4. Set a sufficiently long time span for the integration.
# The solver will stop early if the event is found.
t_span = [0, 20]

# 5. Call the ODE solver.
sol = solve_ivp(
    fun=ode_system,
    t_span=t_span,
    y0=y0,
    events=reach_target_b,
    dense_output=True  # Recommended for accurate event location.
)

# 6. Output the results.
# The final equation we aim to solve is b(t) = 0.5, where t is the unknown.
# The problem gives the numbers for the system: A=1, k=5, a(0)=0.1, b(0)=2.
# The code below will print these numbers and the final computed value for t.
print(f"Solving the system of ODEs with parameters A = {A} and k = {k}.")
print(f"The initial conditions are a(0) = {a0} and b(0) = {b0}.")
print(f"The target is to find the time t that satisfies the equation: b(t) = {b_target}.")
print("-" * 30)

# Check if the event was found and print the result.
if sol.status == 1 and sol.t_events[0].size > 0:
    t_event = sol.t_events[0][0]
    print(f"The time t at which b(t) = {b_target} was found to be: {t_event:.4f}")
else:
    print(f"A time t where b(t) = {b_target} was not found in the interval {t_span}.")
