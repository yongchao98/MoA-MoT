import numpy as np
from scipy.integrate import solve_ivp
import sys

# Suppress warnings for cleaner output, not essential for the logic
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def solve_and_find_time():
    """
    This function sets up and solves the ODE system to find the time t
    when b(t) = 0.5.
    """
    # Step 1: Define the system of ODEs
    def ode_system(t, y, A, k):
        """
        Defines the system of differential equations:
        a'(t) = -0.5*a^2 - A*b^2 + k*(b - 1)
        b'(t) = -a*b
        
        y is a list or array where y[0] = a(t) and y[1] = b(t).
        """
        a, b = y
        dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        dbdt = -a * b
        return [dadt, dbdt]

    # Step 2: Set parameters and initial conditions
    A = 1.0
    k = 5.0
    a0 = 0.1
    b0 = 2.0
    initial_conditions = [a0, b0]

    # Print the problem statement with specified values, including all numbers
    # in the final equations as requested.
    print("The system of differential equations is:")
    print(f"a'(t) = -0.5 * a(t)^2 - {A} * b(t)^2 + {k} * (b(t) - 1)")
    print(f"b'(t) = -a(t) * b(t)")
    print("\nGiven initial conditions:")
    print(f"a(0) = {a0}")
    print(f"b(0) = {b0}")

    # Step 3: Define the event to find when b(t) = 0.5
    target_b_value = 0.5

    def event_b_reaches_target(t, y, A, k):
        """Event function: returns zero when b(t) equals the target value."""
        return y[1] - target_b_value

    # Make the integration stop when the event occurs.
    event_b_reaches_target.terminal = True
    # The event should trigger when b(t) is decreasing.
    event_b_reaches_target.direction = -1

    # Set the time span for integration. A large interval ensures the event is found.
    t_span = [0, 50]

    # Step 4: Solve the ODE system numerically
    solution = solve_ivp(
        ode_system,
        t_span,
        initial_conditions,
        args=(A, k),
        events=event_b_reaches_target,
        dense_output=True 
    )

    # Step 5: Output the result
    print(f"\nSolving for the time t at which b(t) = {target_b_value}...")
    if solution.status == 1 and len(solution.t_events[0]) > 0:
        # A status of 1 indicates the integration was terminated by an event.
        time_of_event = solution.t_events[0][0]
        print(f"\nThe time t at which b(t) = {target_b_value} is approximately: {time_of_event:.4f}")
    else:
        # This case handles if the event is not found in the time span.
        print(f"\nThe condition b(t) = {target_b_value} was not met in the time interval {t_span}.")
        print("This could mean that no such t exists, or a larger time span is needed.")

# Execute the main function
solve_and_find_time()