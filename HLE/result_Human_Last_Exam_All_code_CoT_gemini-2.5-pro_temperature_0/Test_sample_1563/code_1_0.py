import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_time():
    """
    Solves the system of ODEs and finds the time t when b(t) = 0.5.
    """
    # Parameters and initial conditions from the problem
    k = 5
    A = 1
    a0 = 0.1
    b0 = 2
    target_b = 0.5

    # Print the final equations with the numbers substituted in, as requested.
    print("System of differential equations to solve:")
    print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
    print(f"b'(t) = -a(t)*b(t)")
    print(f"\nInitial conditions: a(0) = {a0}, b(0) = {b0}")
    print(f"Goal: Find the time t when b(t) = {target_b}\n")

    # Define the system of ODEs for the solver
    # y is a vector [a(t), b(t)]
    def ode_system(t, y):
        a, b = y
        dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        dbdt = -a * b
        return [dadt, dbdt]

    # Define the event function to find when b(t) = 0.5
    # The solver looks for a root of this function (i.e., when it equals zero)
    def event_b_reaches_target(t, y):
        return y[1] - target_b
    
    # Stop the integration when the event is found
    event_b_reaches_target.terminal = True
    # The event should trigger when b(t) is decreasing
    event_b_reaches_target.direction = -1

    # Set the time span for integration. We choose a large enough interval
    # to ensure the event is found if it exists.
    t_span = [0, 20]
    
    # Initial conditions vector
    y0 = [a0, b0]

    # Call the ODE solver
    solution = solve_ivp(
        ode_system,
        t_span,
        y0,
        events=event_b_reaches_target,
        dense_output=True
    )

    # Extract and print the result
    if solution.t_events[0].size > 0:
        event_time = solution.t_events[0][0]
        print(f"The condition b(t) = {target_b} is met at t = {event_time:.4f}")
    else:
        print(f"The condition b(t) = {target_b} was not met in the time interval {t_span}.")

# Execute the function
solve_and_find_time()