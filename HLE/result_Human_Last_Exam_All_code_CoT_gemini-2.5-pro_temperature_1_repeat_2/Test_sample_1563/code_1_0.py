import numpy as np
from scipy.integrate import solve_ivp

def solve_ode_system():
    """
    This function solves the given system of ODEs to find the time t
    at which b(t) = 0.5.
    """
    # 1. Define the system of differential equations
    def ode_system(t, y, A, k):
        """
        Defines the system of ODEs.
        y[0] represents a(t)
        y[1] represents b(t)
        """
        a, b = y
        da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        db_dt = -a * b
        return [da_dt, db_dt]

    # 2. Set the parameters and initial conditions
    A = 1
    k = 5
    a0 = 0.1
    b0 = 2
    y0 = [a0, b0]
    target_b_value = 0.5

    # 3. Define the event function to find when b(t) = 0.5
    def event_b_reaches_target(t, y, A, k):
        """
        Event function that triggers when b(t) equals the target value.
        The solver finds the root of this function.
        """
        return y[1] - target_b_value

    # Stop the integration when the event occurs
    event_b_reaches_target.terminal = True
    # Trigger the event only when b(t) is decreasing
    event_b_reaches_target.direction = -1

    # 4. Set a time span for the integration.
    # It should be long enough for the condition to be met.
    t_span = [0, 50]

    # 5. Solve the ODE system
    solution = solve_ivp(
        fun=lambda t, y: ode_system(t, y, A, k),
        t_span=t_span,
        y0=y0,
        events=lambda t, y: event_b_reaches_target(t, y, A, k),
        dense_output=True  # Needed to evaluate the solution at the event time
    )

    # 6. Output the results
    print("Solving the system of differential equations:")
    # The final code outputs each number in the final equation as requested.
    print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t) - 1)")
    print(f"b'(t) = -a(t)*b(t)")
    print(f"with initial conditions a(0) = {a0}, b(0) = {b0}.")
    print("\nFinding the time t such that b(t) = 0.5...\n")


    # Check if the event was found
    if solution.t_events[0].size > 0:
        t_at_event = solution.t_events[0][0]
        state_at_event = solution.sol(t_at_event)
        a_at_event = state_at_event[0]
        b_at_event = state_at_event[1]
        
        print(f"Result:")
        print(f"The time t at which b(t) = {target_b_value} is approximately: {t_at_event:.4f}")
        print(f"The state of the system at this time is:")
        print(f"a({t_at_event:.4f}) = {a_at_event:.4f}")
        print(f"b({t_at_event:.4f}) = {b_at_event:.4f}")
    else:
        print(f"The condition b(t) = {target_b_value} was not met in the time interval t = [0, {t_span[1]}].")

if __name__ == "__main__":
    solve_ode_system()