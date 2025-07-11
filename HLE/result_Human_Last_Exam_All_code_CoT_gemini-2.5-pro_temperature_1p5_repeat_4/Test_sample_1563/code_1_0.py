import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_time():
    """
    Solves the system of ODEs and finds the time t when b(t) reaches a target value.
    """
    # --- Problem Setup ---
    # The system of differential equations is:
    # a'(t) = -1/2 * a^2(t) - A * b^2(t) + k * (b(t) - 1)
    # b'(t) = -a(t) * b(t)

    # --- Parameters and Initial Conditions ---
    k = 5.0
    A = 1.0
    a0 = 0.1
    b0 = 2.0
    target_b = 0.5
    
    initial_conditions = [a0, b0]

    # --- Model Definition ---
    # We define the model function for the solver.
    # The state vector y is [a(t), b(t)].
    def differential_equations_model(t, y):
        """
        Defines the system of ODEs.
        t: time
        y: state vector [a, b]
        Returns: [da/dt, db/dt]
        """
        a, b = y
        dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        dbdt = -a * b
        return [dadt, dbdt]

    # --- Event Definition ---
    # We need to find the time 't' where b(t) equals target_b.
    # We use the 'events' feature of solve_ivp to find this time precisely.
    def event_b_reaches_target(t, y):
        """
        Event function that triggers when b(t) - target_b = 0.
        y[1] corresponds to b(t).
        """
        return y[1] - target_b

    # This option stops the integration when the event is found.
    event_b_reaches_target.terminal = True
    # We expect b(t) to decrease, so we trigger the event when the function
    # value is decreasing and crosses zero.
    event_b_reaches_target.direction = -1

    # --- Solving the ODE ---
    # Set a long enough time span for the integration to ensure the event can be found.
    t_span = [0, 50]

    # Call the numerical solver.
    solution = solve_ivp(
        fun=differential_equations_model,
        t_span=t_span,
        y0=initial_conditions,
        events=event_b_reaches_target,
        dense_output=True 
    )

    # --- Output the Result ---
    print("Solving the system of differential equations:")
    print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
    print(f"b'(t) = -a(t)*b(t)")
    print(f"with initial conditions a(0) = {a0}, b(0) = {b0}.")
    print(f"We want to find the time t at which b(t) = {target_b}.")
    print("-" * 40)

    # Check if the solver found the event and print the corresponding time.
    if solution.t_events[0].size > 0:
        time_at_event = solution.t_events[0][0]
        print(f"The condition b(t) = {target_b} was met at time t = {time_at_event:.4f}")
    else:
        # This branch handles the case where the integration finishes without the event happening.
        print(f"The condition b(t) = {target_b} was not met in the time interval from t={t_span[0]} to t={t_span[1]}.")

if __name__ == '__main__':
    solve_and_find_time()