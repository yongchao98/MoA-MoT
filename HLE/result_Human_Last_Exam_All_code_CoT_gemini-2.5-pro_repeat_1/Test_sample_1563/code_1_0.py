import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_time():
    """
    Solves the system of ODEs to find the time t when b(t) = 0.5.
    """
    # Define the system of differential equations
    # y[0] = a(t), y[1] = b(t)
    def model(t, y):
        a, b = y
        # Constants from the problem statement
        k = 5
        A = 1
        
        # The equations
        # a'(t) = -1/2 * a^2(t) - A * b^2(t) + k * (b(t) - 1)
        # b'(t) = -a(t) * b(t)
        dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        dbdt = -a * b
        return [dadt, dbdt]

    # Set the initial conditions: (a(0), b(0)) = (0.1, 2)
    initial_conditions = [0.1, 2]

    # Set the time span for the integration.
    # Choose a time long enough for the event to occur.
    t_span = (0, 20)

    # Define the event: stop when b(t) = 0.5
    # The event function should be zero at the event time.
    def event_b_reaches_half(t, y):
        return y[1] - 0.5
    
    # The integration should stop when the event is found.
    event_b_reaches_half.terminal = True
    # From the initial conditions, b(t) is decreasing, so the direction is negative.
    event_b_reaches_half.direction = -1

    # Solve the IVP
    solution = solve_ivp(
        fun=model,
        t_span=t_span,
        y0=initial_conditions,
        events=event_b_reaches_half,
        dense_output=True
    )

    # Check if the event was found and print the result
    if solution.t_events and solution.t_events[0].size > 0:
        event_time = solution.t_events[0][0]
        # The final equation is b(t) = 0.5, which is satisfied at the calculated time.
        # We print out the number in this final state.
        b_final = 0.5
        print(f"The time t at which b(t) = {b_final} is: {event_time}")
    else:
        print("The event b(t) = 0.5 was not found in the specified time interval.")

# Run the solver
solve_and_find_time()