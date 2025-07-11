import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_time():
    """
    This function solves the given system of ODEs and finds the time t
    at which b(t) = 0.5.
    """
    # Define the system of differential equations
    def system_of_odes(t, y):
        """
        Defines the system of ODEs.
        y[0] = a(t)
        y[1] = b(t)
        """
        a, b = y
        k = 5.0
        A = 1.0
        
        # a'(t) = -1/2 * a^2(t) - A * b^2(t) + k * (b(t) - 1)
        dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        
        # b'(t) = -a(t) * b(t)
        dbdt = -a * b
        
        return [dadt, dbdt]

    # Set the initial conditions
    # (a(0), b(0)) = (0.1, 2)
    initial_conditions = [0.1, 2.0]

    # Set the time span for integration. A sufficiently large end time is chosen.
    t_span = [0, 20]

    # Define the event function to find when b(t) = 0.5
    def event_b_equals_0_5(t, y):
        """
        Event function that triggers when b(t) = 0.5.
        We are looking for the root of y[1] - 0.5 = 0.
        """
        return y[1] - 0.5

    # Configure the event to stop the integration when it occurs.
    # Since b(0)=2, b(t) is decreasing, so we look for a negative-going crossing.
    event_b_equals_0_5.terminal = True
    event_b_equals_0_5.direction = -1

    # Solve the ODE system using the 'RK45' method, which is a good default.
    # dense_output=True allows for accurate evaluation of the solution between time steps.
    solution = solve_ivp(
        system_of_odes,
        t_span,
        initial_conditions,
        events=event_b_equals_0_5,
        dense_output=True
    )

    # Check if the event was triggered and print the result.
    if solution.status == 1: # status=1 means an event terminated the integration
        t_event = solution.t_events[0][0]
        # To fulfill the request "output each number in the final equation",
        # we print the equation b(t) = 0.5 with the found value of t.
        # We also retrieve the value of b(t) at the event time to confirm it is 0.5.
        y_event = solution.sol(t_event)
        b_event = y_event[1]
        print(f"The time t at which b(t) reaches 0.5 was found.")
        print(f"Final Equation: b({t_event:.4f}) = {b_event:.4f}")
    else:
        print("The condition b(t) = 0.5 was not met in the specified time interval.")

# Execute the function to find the answer
solve_and_find_time()