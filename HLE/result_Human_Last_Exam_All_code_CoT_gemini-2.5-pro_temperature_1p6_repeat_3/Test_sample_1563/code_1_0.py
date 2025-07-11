import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_time():
    """
    Solves the given system of ODEs and finds the time t when b(t) = 0.5.
    """
    # --- 1. Define the ODE system ---
    # Parameters are defined inside the function to keep it self-contained.
    A = 1.0
    k = 5.0

    def ode_system(t, y):
        """
        Defines the system of differential equations.
        y[0] = a(t)
        y[1] = b(t)
        """
        a, b = y
        da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        db_dt = -a * b
        return [da_dt, db_dt]

    # --- 2. Set up the problem ---
    # Initial conditions
    a0 = 0.1
    b0 = 2.0
    y0 = [a0, b0]

    # Time interval for integration. A sufficiently large interval.
    t_span = [0, 20]

    # --- 3. Define the event to find t when b(t)=0.5 ---
    def event_b_equals_0_5(t, y):
        """Event function, returns 0 when b(t) = 0.5"""
        return y[1] - 0.5

    # Stop integration when this event occurs
    event_b_equals_0_5.terminal = True
    # The event should be triggered when b(t) is decreasing towards 0.5
    event_b_equals_0_5.direction = -1

    # --- 4. Solve the ODE system ---
    sol = solve_ivp(
        fun=ode_system,
        t_span=t_span,
        y0=y0,
        events=event_b_equals_0_5,
        dense_output=True # Needed for accurate event location
    )

    # --- 5. Output the results ---
    print("The system of differential equations with k=5 and A=1 is:")
    print(f"a'(t) = -0.5*a(t)^2 - {A:.1f}*b(t)^2 + {k:.1f}*(b(t)-1)")
    print(f"b'(t) = -a(t)*b(t)")
    print(f"\nWith initial conditions: a(0) = {a0}, b(0) = {b0}")

    if sol.status == 1 and len(sol.t_events[0]) > 0:
        # The event was found
        t_at_event = sol.t_events[0][0]
        print(f"\nThe time t at which b(t) equals 0.5 is approximately {t_at_event:.4f}.")
    else:
        # The event was not found in the integration interval
        print("\nThe condition b(t) = 0.5 was not met in the specified time interval.")

# Run the solver
solve_and_find_time()