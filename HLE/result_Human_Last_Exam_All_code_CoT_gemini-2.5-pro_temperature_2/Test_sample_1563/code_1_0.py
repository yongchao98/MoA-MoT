import numpy as np
from scipy.integrate import solve_ivp

def solve_and_find_time():
    """
    Solves the given system of ODEs to find the time t when b(t) = 0.5.
    """
    # 1. Define the parameters from the problem description
    k = 5.0
    A = 1.0

    # 2. Define the system of differential equations
    # y[0] represents a(t) and y[1] represents b(t)
    def model(t, y):
        a, b = y
        # a'(t) = -1/2 * a^2(t) - A * b^2(t) + k * (b(t)-1)
        dadt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        # b'(t) = -a(t) * b(t)
        dbdt = -a * b
        return [dadt, dbdt]

    # 3. Set the initial conditions
    a0 = 0.1
    b0 = 2.0
    initial_conditions = [a0, b0]

    # 4. Define the time span for the integration.
    # We choose a large enough end time to ensure the event is found.
    t_span = (0, 20)

    # 5. Define the event function to find t when b(t) = 0.5
    # The solver looks for a zero crossing of this function.
    def reach_target_b(t, y):
        return y[1] - 0.5

    # Make the integration terminate when the event occurs
    reach_target_b.terminal = True
    # The event should trigger only when b(t) is decreasing
    reach_target_b.direction = -1

    # 6. Call the ODE solver
    solution = solve_ivp(
        fun=model,
        t_span=t_span,
        y0=initial_conditions,
        events=reach_target_b,
        dense_output=True  # for accurate event location
    )

    # 7. Output the results
    print("Solving the system of differential equations:")
    print(f"a'(t) = -{0.5}*a^2(t) - {A}*b^2(t) + {k}*(b(t)-1)")
    print(f"b'(t) = -a(t)*b(t)")
    print(f"\nGiven initial conditions: a(0) = {a0}, b(0) = {b0}")
    print(f"Target condition: b(t) = 0.5")

    # Check if the event was found and print the time
    if solution.status == 1 and solution.t_events[0].size > 0:
        event_time = solution.t_events[0][0]
        print(f"\nEstimated time t at which b(t) = 0.5 is: {event_time:.4f}")
    else:
        print("\nThe target condition was not met in the specified time interval.")

if __name__ == "__main__":
    solve_and_find_time()