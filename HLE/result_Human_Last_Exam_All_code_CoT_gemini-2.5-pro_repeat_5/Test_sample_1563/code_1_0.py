import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_time():
    """
    Solves the system of ODEs and finds the time t when b(t) = 0.5.
    """
    # 1. Define the system of differential equations
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

    # 2. Set parameters and initial conditions
    k = 5.0
    A = 1.0
    a0 = 0.1
    b0 = 2.0
    y0 = [a0, b0]

    # 3. Define the event to find when b(t) = 0.5
    # The event function should return 0 at the event.
    def event_b_reaches_0_5(t, y, A, k):
        return y[1] - 0.5

    # The solver should stop when this event is detected.
    event_b_reaches_0_5.terminal = True
    # The event function is decreasing (b(t) goes from 2 towards 0.5).
    event_b_reaches_0_5.direction = -1

    # 4. Set a sufficiently long time span for the integration
    t_span = [0, 20]

    # 5. Solve the ODE system
    sol = solve_ivp(
        fun=model,
        t_span=t_span,
        y0=y0,
        args=(A, k),
        events=event_b_reaches_0_5,
        dense_output=True # Allows for accurate interpolation if needed
    )

    # 6. Output the results
    print("Solving the system of differential equations:")
    print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
    print(f"b'(t) = -a(t)*b(t)")
    print(f"\nWith initial conditions a(0) = {a0}, b(0) = {b0}.")
    
    if sol.t_events[0].size > 0:
        t_event = sol.t_events[0][0]
        b_event = sol.y_events[0][0][1]
        print(f"\nThe condition b(t) = {b_event:.1f} is met at time t = {t_event:.4f}")
    else:
        print(f"\nThe condition b(t) = 0.5 was not reached in the time interval {t_span}.")

# Run the solver
solve_and_find_time()