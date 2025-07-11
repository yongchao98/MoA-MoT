import numpy as np
from scipy.integrate import solve_ivp

# Step 1: Define the system of differential equations
def model(t, y, A, k):
    """
    Defines the system of ODEs.
    y[0] = a(t)
    y[1] = b(t)
    """
    a, b = y
    dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
    dbdt = -a * b
    return [dadt, dbdt]

# Step 2: Define the constants and initial conditions
A = 1
k = 5
# Initial conditions: a(0) = 0.1, b(0) = 2
y0 = [0.1, 2]

# Step 3: Define the event function to find when b(t) = 0.5
def event_b_equals_0_5(t, y, A, k):
    """
    Event function that triggers when b(t) - 0.5 = 0.
    """
    return y[1] - 0.5

# We expect b(t) to decrease to 0.5, so we specify the direction of the zero-crossing.
event_b_equals_0_5.direction = -1
# We want the integration to stop when the event occurs.
event_b_equals_0_5.terminal = True

# Step 4: Set the time span for the integration
# Choose an end time large enough for the event to occur.
t_span = [0, 20]

# Step 5: Solve the ODE system using solve_ivp
sol = solve_ivp(
    model,
    t_span,
    y0,
    args=(A, k),
    events=event_b_equals_0_5,
    dense_output=True
)

# Step 6: Extract and print the results
print(f"The system of differential equations is:")
print(f"a'(t) = -0.5 * a(t)^2 - {A} * b(t)^2 + {k} * (b(t) - 1)")
print(f"b'(t) = -a(t) * b(t)")
print(f"\nGiven initial conditions: a(0) = {y0[0]}, b(0) = {y0[1]}")
print(f"We are looking for the time t at which b(t) = 0.5.")

# Check if the event was found
if sol.status == 1 and len(sol.t_events[0]) > 0:
    # The time at which the event occurred
    t_event = sol.t_events[0][0]
    # The state of the system at that time
    a_event, b_event = sol.y_events[0][0]

    print(f"\nNumeric integration finds that b(t) reaches 0.5 at time t = {t_event:.4f}.")
    print(f"The state at this time is a(t) = {a_event:.4f} and b(t) = {b_event:.4f}.")

    print("\nAt this point in time, the values of the derivatives ('the final equation') are:")
    # Calculate the derivatives at the event time
    dadt_event, dbdt_event = model(t_event, [a_event, b_event], A, k)

    # Print the detailed calculation for a'(t)
    print(f"a'({t_event:.2f}) = -0.5 * ({a_event:.4f})^2 - {A} * ({b_event:.4f})^2 + {k} * ({b_event:.4f} - 1) = {dadt_event:.4f}")

    # Print the detailed calculation for b'(t)
    print(f"b'({t_event:.2f}) = -({a_event:.4f}) * ({b_event:.4f}) = {dbdt_event:.4f}")
else:
    print(f"\nThe condition b(t) = 0.5 was not met in the time interval from t={t_span[0]} to t={t_span[1]}.")
    print("This implies no such t exists or it is larger than the integration time.")
