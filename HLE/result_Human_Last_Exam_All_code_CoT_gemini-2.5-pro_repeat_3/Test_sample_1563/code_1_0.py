import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_time():
    """
    Solves the system of ODEs and finds the time t when b(t) = 0.5.
    """
    # Define the system of differential equations
    # y[0] = a(t), y[1] = b(t)
    def ode_system(t, y, A, k):
        a, b = y
        dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        dbdt = -a * b
        return [dadt, dbdt]

    # Given parameters and initial conditions
    A = 1
    k = 5
    a0 = 0.1
    b0 = 2
    y0 = [a0, b0]

    # Time span for the integration. We choose a large enough interval
    # to ensure the event is found.
    t_span = [0, 20]

    # Define the event function to find t when b(t) = 0.5
    # The solver looks for a zero of this function.
    def event_b_reaches_0_5(t, y):
        return y[1] - 0.5

    # Make the event terminal: the integration will stop when the event is found.
    event_b_reaches_0_5.terminal = True
    # The event function is decreasing (b(t) goes from 2 towards 0.5).
    event_b_reaches_0_5.direction = -1

    # Solve the ODE system
    sol = solve_ivp(
        fun=lambda t, y: ode_system(t, y, A, k),
        t_span=t_span,
        y0=y0,
        events=event_b_reaches_0_5,
        dense_output=True
    )

    # Print the final equations with the given parameter values
    print("The system of differential equations being solved is:")
    # The prompt asks to "output each number in the final equation"
    print(f"a'(t) = -1/2 * a^2(t) - {A} * b^2(t) + {k} * (b(t) - 1)")
    print(f"b'(t) = -a(t) * b(t)")
    print(f"\nInitial conditions: a(0) = {a0}, b(0) = {b0}")
    print("\nThe goal is to find the time t at which b(t) = 0.5.\n")


    # Check if the event was found and print the result
    if sol.status == 1 and sol.t_events[0].size > 0:
        t_event = sol.t_events[0][0]
        print(f"The event b(t) = 0.5 occurred at time t = {t_event:.4f}")
    else:
        print("The event b(t) = 0.5 was not found in the integration interval.")

if __name__ == '__main__':
    solve_and_find_time()